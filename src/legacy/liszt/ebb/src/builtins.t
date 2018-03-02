-- The MIT License (MIT)
--
-- Copyright (c) 2015 Stanford University.
-- All rights reserved.
--
-- Permission is hereby granted, free of charge, to any person obtaining a
-- copy of this software and associated documentation files (the "Software"),
-- to deal in the Software without restriction, including without limitation
-- the rights to use, copy, modify, merge, publish, distribute, sublicense,
-- and/or sell copies of the Software, and to permit persons to whom the
-- Software is furnished to do so, subject to the following conditions:
--
-- The above copyright notice and this permission notice shall be included
-- in all copies or substantial portions of the Software.
--
-- THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
-- IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
-- FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
-- AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
-- LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
-- FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
-- DEALINGS IN THE SOFTWARE.
local B = {}
package.loaded["ebb.src.builtins"] = B

local Pre = require "ebb.src.prelude"
local T   = require "ebb.src.types"
local C   = require "ebb.src.c"
local AST = require "ebb.src.ast"
local R   = require "ebb.src.relations"


local errorT    = T.error
local floatT    = T.float
local doubleT   = T.double
local intT      = T.int
local uint64T   = T.uint64
local boolT     = T.bool

local keyT      = T.key
local vectorT   = T.vector
--local matrixT   = T.matrix

local CPU       = Pre.CPU


---------------------------------------------
--[[ Builtin functions                   ]]--
---------------------------------------------
local Builtin = {}
Builtin.__index = Builtin
Builtin.__call = function(self,...)
    return self.luafunc(...)
end

function Builtin.new(luafunc)
    local check = function(ast, ctxt)
        error('unimplemented builtin function typechecking')
    end
    if not luafunc then
        luafunc = function() error("Cannot call from Lua code") end
    end
    return setmetatable({check=check, luafunc = luafunc},
                        Builtin)
end
function B.is_builtin(f)
    return getmetatable(f) == Builtin
end

local function id_checks(fname, ast, ctxt, args)
    if #args ~= 1 then
        ctxt:error(ast, fname.." expects exactly 1 argument (instead got " ..
                        tostring(#args) .. ")")
        return false
    end

    if not args[1].node_type:isscalarkey() then
        ctxt:error(ast, "expected a relational key as the argument for "..
                        fname.."()")
        return false
    end
    return true
end

B.id = Builtin.new()
function B.id.check(ast, ctxt)
    local args = ast.params

    if not id_checks('id', ast, ctxt, args) then return errorT end
    if args[1].node_type.ndims ~= 1 then
        ctxt:error(ast, "Can only use built-in id() on keys of "..
                        "non-grid relations; "..
                        "try using xid(), yid() or zid() instead.")
        return errorT
    end
    if args[1].node_type.relation:isCoupled() then
        ctxt:error(ast, "Can't access the id() of a coupled relation.")
        return errorT
    end


    return uint64T
end

B.xid = Builtin.new()
function B.xid.check(ast, ctxt)
    local args = ast.params

    if not id_checks('xid', ast, ctxt, args) then return errorT end
    if args[1].node_type.ndims == 1 then
        ctxt:error(ast, "Can only use built-in xid() on keys of "..
                        "grid relations; try using id() instead.")
        return errorT
    end

    return uint64T
end

B.yid = Builtin.new()
function B.yid.check(ast, ctxt)
    local args = ast.params

    if not id_checks('yid', ast, ctxt, args) then return errorT end
    if args[1].node_type.ndims == 1 then
        ctxt:error(ast, "Can only use built-in yid() on keys of "..
                        "grid relations; try using id() instead.")
        return errorT
    end

    return uint64T
end

B.zid = Builtin.new()
function B.zid.check(ast, ctxt)
    local args = ast.params

    if not id_checks('zid', ast, ctxt, args) then return errorT end
    if args[1].node_type.ndims == 1 then
        ctxt:error(ast, "Can only use built-in zid() on keys of "..
                        "grid relations; try using id() instead.")
        return errorT
    end
    if args[1].node_type.ndims < 3 then
        ctxt:error(ast, "The key argument to zid() refers to a 2d grid, "..
                        "so zid() doesn't make any sense.")
        return errorT
    end

    return uint64T
end


B.Affine = Builtin.new()
function B.Affine.check(ast, ctxt)
    local args = ast.params

    if #args ~= 3 then
        ctxt:error(ast,'Affine expects 3 arguments')
        return errorT
    end
    local dst_rel_arg   = args[1]
    local matrix        = args[2]
    local key_arg       = args[3]
    local ret_type      = nil

    -- check that the first and last arg are actually relations
    if not dst_rel_arg.node_type:isinternal() or
       not R.is_relation(dst_rel_arg.node_type.value)
    then
        ctxt:error(ast[1], "Affine expects a relation as the 1st argument")
        return errorT
    end
    if not key_arg.node_type:isscalarkey() then
        ctxt:error(ast[3], "Affine expects a key as the 3rd argument")
        return errorT
    end

    -- get the source and destination relations and check that they're grids
    local dst_rel = dst_rel_arg.node_type.value
    local src_rel = key_arg.node_type.relation
    if not dst_rel:isGrid() then
        ctxt:error(ast[1],
            "Affine expects a grid relation as the 1st argument")
        return errorT
    end
    if not src_rel:isGrid() then
        ctxt:error(ast[3], "Affine expects a grid key as the 3rd argument")
        return errorT
    end

    -- get dimensions out
    local dst_dims = dst_rel:NumDims()
    local src_dims = src_rel:NumDims()

    -- now check the matrix argument type
    if not matrix.node_type:ismatrix() or
       matrix.node_type.Nrow ~= dst_dims or
       matrix.node_type.Ncol ~= src_dims + 1
    then
        ctxt:error(ast[2], "Affine expects a matrix as the 2nd argument "..
            "with matching dimensions (needs to be "..
            tostring(dst_dims).."-by-"..tostring(src_dims + 1))
        return errorT
    end
    --if not matrix.node_type:isintegral() then
    --    ctxt:error(ast[2], "Affine expects a matrix of integral values")
    --    return errorT
    --end
    if not matrix.node_type:isnumeric() then
        ctxt:error(ast[2], "Affine expects a matrix of numeric values")
        return errorT
    end
    -- Opting for a weak literal test instead of a const-ness test for now
    if not matrix:is(AST.MatrixLiteral) or
       not matrix:isLiteral()
    then
        ctxt:error(args[2], "Compiler could not verify that "..
            "the matrix argument (2nd) to Affine is a literal")
        return errorT
    end

    -- cache the matrix literal value for easy access
    local Nrow,Ncol = matrix.node_type.Nrow, matrix.node_type.Ncol
    matrix.matvals   = {}
    for i=1,Nrow do
        matrix.matvals[i] = {}
        for j=1,Ncol do
            local entry_ast = matrix.elems[(i-1)*Ncol + j]
            local val = assert(entry_ast:is(AST.Number) and entry_ast.value,
                               'INTERNAL: found non-literal matrix entry')
            matrix.matvals[i][j] = val
        end
    end

    return keyT(dst_rel)
end
local terra full_mod(val : int64, modulus : int64) : uint64
    return ((val % modulus) + modulus) % modulus
end

B.UNSAFE_ROW = Builtin.new()
function B.UNSAFE_ROW.check(ast, ctxt)
    local args = ast.params
    if #args ~= 2 then
        ctxt:error(ast, "UNSAFE_ROW expects exactly 2 arguments "..
                        "(instead got " .. #args .. ")")
        return errorT
    end

    local addr_type = args[1].node_type
    local rel_type = args[2].node_type
    if not rel_type:isinternal() or not R.is_relation(rel_type.value) then
        ctxt:error(ast, "UNSAFE_ROW expected a relation as the second arg")
        return errorT
    end
    local rel = rel_type.value
    if rel:isCoupled() then
        ctxt:error(ast, "UNSAFE_ROW can't be used for coupled relations.")
        return errorT
    end
    local ndim = rel:NumDims()
    if ndim == 1 and addr_type ~= uint64T then
        ctxt:error(ast, "UNSAFE_ROW expected a uint64 as the first arg")
        return errorT
    elseif ndim > 1  and addr_type ~= vectorT(uint64T,ndim) then
        ctxt:error(ast, "UNSAFE_ROW expected a vector of "..ndim..
                        " uint64 values")
        return errorT
    end

    -- actual typing
    return keyT(rel)
end


B.assert = Builtin.new(assert)
function B.assert.check(ast, ctxt)
    local args = ast.params
    if #args ~= 1 then
        ctxt:error(ast, "assert expects exactly 1 argument (instead got " .. #args .. ")")
        return
    end

    local test = args[1]
    local test_type = test.node_type
    if test_type:isvector() then test_type = test_type:basetype() end
    if test_type ~= errorT and test_type ~= boolT then
        ctxt:error(ast, "expected a boolean or vector of booleans as the test for assert statement")
    end
end

local terra ebbAssert(test : bool, file : rawstring, line : int)
    if not test then
        C.fprintf(C.stderr, "%s:%d: assertion failed!\n", file, line)
        C.exit(1)
    end
end

B.print = Builtin.new(print)
function B.print.check(ast, ctxt)
    local args = ast.params

    for i,output in ipairs(args) do
        local outtype = output.node_type
        if outtype ~= errorT and
           not outtype:isvalue() and not outtype:iskey()
        then
            ctxt:error(ast, "only numbers, bools, vectors, matrices and keys can be printed")
        end
    end
end

local function printSingle (bt, exp, elemQuotes)
    if bt == floatT or bt == doubleT then
        table.insert(elemQuotes, `[double]([exp]))
        return "%f"
    elseif bt == intT then
        table.insert(elemQuotes, exp)
        return "%d"
    elseif bt == uint64T then
        table.insert(elemQuotes, exp)
        return "%lu"
    elseif bt == boolT then
        table.insert(elemQuotes, `terralib.select([exp], "true", "false"))
        return "%s"
    else
        error('Unrecognized type in print: ' .. bt:toString() .. ' ' .. tostring(bt:terratype()))
    end
end

B.rand = Builtin.new()
function B.rand.check(ast, ctxt)
    local args = ast.params
    if #args ~= 0 then
        ctxt:error(ast, "rand expects 0 arguments")
        return errorT
    else
        return doubleT
    end
end

B.dot = Builtin.new()
function B.dot.check(ast, ctxt)
    local args = ast.params
    if #args ~= 2 then
        ctxt:error(ast, "dot product expects exactly 2 arguments "..
                        "(instead got " .. #args .. ")")
        return errorT
    end

    local lt1 = args[1].node_type
    local lt2 = args[2].node_type

    local numvec_err = 'arguments to dot product must be numeric vectors'
    local veclen_err = 'vectors in dot product must have equal dimensions'
    if not lt1:isvector()  or not lt2:isvector() or
       not lt1:isnumeric() or not lt2:isnumeric()
    then
        ctxt:error(ast, numvec_err)
    elseif lt1.N ~= lt2.N then
        ctxt:error(ast, veclen_err)
    else
        return T.type_join(lt1:basetype(), lt2:basetype())
    end

    return errorT
end

B.times = Builtin.new()
function B.times.check(ast, ctxt)
    local args = ast.params
    if #args ~= 2 then
        ctxt:error(ast, "element-wise multiplication expects exactly 2 arguments "..
                        "(instead got " .. #args .. ")")
        return errorT
    end

    local lt1 = args[1].node_type
    local lt2 = args[2].node_type

    local numvec_err = 'arguments to element-wise multiplication must be numeric vectors'
    local veclen_err = 'vectors in element-wise multiplication must have equal dimensions'
    if not lt1:isvector()  or not lt2:isvector() or
       not lt1:isnumeric() or not lt2:isnumeric()
    then
        ctxt:error(ast, numvec_err)
    elseif lt1.N ~= lt2.N then
        ctxt:error(ast, veclen_err)
    else
        return T.type_join(lt1, lt2)
    end

    return errorT
end

B.cross = Builtin.new()
function B.cross.check(ast, ctxt)
    local args = ast.params

    if #args ~= 2 then
        ctxt:error(ast, "cross product expects exactly 2 arguments "..
                        "(instead got " .. #args .. ")")
        return errorT
    end

    local lt1 = args[1].node_type
    local lt2 = args[2].node_type

    local numvec_err = 'arguments to cross product must be numeric vectors'
    local veclen_err = 'vectors in cross product must be 3 dimensional'
    if not lt1:isvector()  or not lt2:isvector() or
       not lt1:isnumeric() or not lt2:isnumeric()
    then
        ctxt:error(ast, numvec_err)
    elseif lt1.N ~= 3 or lt2.N ~= 3 then
        ctxt:error(ast, veclen_err)
    else
        return T.type_join(lt1, lt2)
    end

    return errorT
end

B.length = Builtin.new()
function B.length.check(ast, ctxt)
    local args = ast.params
    if #args ~= 1 then
        ctxt:error(ast, "length expects exactly 1 argument (instead got " .. #args .. ")")
        return errorT
    end
    local lt = args[1].node_type
    if not lt:isvector() then
        ctxt:error(args[1], "argument to length must be a vector")
        return errorT
    end
    if not lt:basetype():isnumeric() then
        ctxt:error(args[1], "length expects vectors of numeric type")
    end
    if lt:basetype() == floatT then return floatT
                                else return doubleT end
end

function Builtin.newDoubleFunction(name)
    local cpu_fn = C[name]
    local lua_fn = function (arg) return cpu_fn(arg) end

    local b = Builtin.new(lua_fn)

    function b.check (ast, ctxt)
        local args = ast.params
        if #args ~= 1 then
            ctxt:error(ast, name.." expects exactly 1 argument "..
                            "(instead got ".. #args ..")")
            return errorT
        end
        local lt = args[1].node_type
        if not lt:isnumeric() then
            ctxt:error(args[1], "argument to "..name.." must be numeric")
        end
        if lt:isvector() then
            ctxt:error(args[1], "argument to "..name.." must be a scalar")
            return errorT
        end
        return doubleT
    end

    function b.codegen (ast, ctxt)
        local exp = ast.params[1]:codegen(ctxt)
          return `cpu_fn([exp])
    end
    return b
end

B.cos   = Builtin.newDoubleFunction('cos')
B.acos  = Builtin.newDoubleFunction('acos')
B.sin   = Builtin.newDoubleFunction('sin')
B.asin  = Builtin.newDoubleFunction('asin')
B.tan   = Builtin.newDoubleFunction('tan')
B.atan  = Builtin.newDoubleFunction('atan')
B.sqrt  = Builtin.newDoubleFunction('sqrt')
B.cbrt  = Builtin.newDoubleFunction('cbrt')
B.floor = Builtin.newDoubleFunction('floor')
B.ceil  = Builtin.newDoubleFunction('ceil')
B.fabs  = Builtin.newDoubleFunction('fabs')
B.log   = Builtin.newDoubleFunction('log')

B.fmin  = Builtin.new()
B.fmax  = Builtin.new()
local function fminmax_check(ast, ctxt, name)
    if #ast.params ~= 2 then
        ctxt:error(ast, name.." expects 2 arguments "..
                        "(instead got ".. #ast.params ..")")
        return errorT
    end
    local lt = ast.params[1].node_type
    local rt = ast.params[1].node_type
    if not lt:isnumeric() or not rt:isscalar() then
        ctxt:error(ast.params[1], "argument to "..name..
                                  " must be a scalar number")
    end
    if not rt:isnumeric() or not rt:isscalar() then
        ctxt:error(ast.params[2], "argument to "..name..
                                  " must be a scalar number")
    end
    return doubleT
end
function B.fmin.check(ast, ctxt) return fminmax_check(ast, ctxt, 'fmin') end
function B.fmax.check(ast, ctxt) return fminmax_check(ast, ctxt, 'fmax') end

B.imin  = Builtin.new()
B.imax  = Builtin.new()
local function iminmax_check(ast, ctxt, name)
    if #ast.params ~= 2 then
        ctxt:error(ast, name.." expects 2 arguments "..
                        "(instead got ".. #ast.params ..")")
        return errorT
    end
    local lt = ast.params[1].node_type
    local rt = ast.params[1].node_type
    if not lt == intT then
        ctxt:error(ast.params[1], "argument to "..name..
                                  " must be an int")
    end
    if not rt == intT then
        ctxt:error(ast.params[2], "argument to "..name..
                                  " must be an int")
    end
    return intT
end
function B.imin.check(ast, ctxt) return iminmax_check(ast, ctxt, 'imin') end
function B.imax.check(ast, ctxt) return iminmax_check(ast, ctxt, 'imax') end

B.pow = Builtin.new(C.pow)
function B.pow.check (ast, ctxt)
    local args = ast.params
    if #args ~= 2 then ctxt:error(ast, "pow expects 2 arguments (instead got " .. #args .. ")")
        return errorT
    end
    for i = 1, #args do
        local lt = args[i].node_type
        if not lt:isnumeric() then
            ctxt:error(args[i], "argument "..i.." to pow must be numeric")
            return errorT
        end
    end
    for i = 1, #args do
        local lt = args[i].node_type
        if not lt:isscalar() then
            ctxt:error(args[i], "argument "..i.." to pow must be a scalar")
            return errorT
        end
    end
    return doubleT
end

B.fmod = Builtin.new(C.fmod)
function B.fmod.check (ast, ctxt)
    local args = ast.params
    if #args ~= 2 then ctxt:error(ast, "fmod expects 2 arguments (instead got " .. #args .. ")")
        return errorT
    end
    for i = 1, #args do
        local lt = args[i].node_type
        if not lt:isnumeric() then
            ctxt:error(args[i], "argument "..i.." to fmod must be numeric")
            return errorT
        end
    end
    for i = 1, #args do
        local lt = args[i].node_type
        if not lt:isscalar() then
            ctxt:error(args[i], "argument "..i.." to fmod must be a scalar")
            return errorT
        end
    end
    return doubleT
end


B.all = Builtin.new()
function B.all.check(ast, ctxt)
    local args = ast.params
    if #args ~= 1 then
        ctxt:error(ast, "all expects exactly 1 argument (instead got " .. #args .. ")")
        return errorT
    end
    local lt = args[1].node_type
    if not lt:isvector() then
        ctxt:error(args[1], "argument to all must be a vector")
        return errorT
    end
    return boolT
end


B.any = Builtin.new()
function B.any.check(ast, ctxt)
    local args = ast.params
    if #args ~= 1 then
        ctxt:error(ast, "any expects exactly 1 argument (instead got " .. #args .. ")")
        return errorT
    end
    local lt = args[1].node_type
    if not lt:isvector() then
        ctxt:error(args[1], "argument to any must be a vector")
        return errorT
    end
    return boolT
end


local function map(fn, list)
    local result = {}
    for i = 1, #list do
        table.insert(result, fn(list[i]))
    end
    return result
end

local function GetTypedSymbol(arg)

    return symbol(arg.node_type:terratype())
end

local function TerraCheck(func)
    return function (ast, ctxt)
        local args = ast.params
        local argsyms = map(GetTypedSymbol, args)
        local rettype = nil
        local try = function()
            local terrafunc = terra([argsyms]) return func([argsyms]) end
            rettype = terrafunc:gettype().returntype
        end
        local success, retval = pcall(try)
        if not success then
            ctxt:error(ast, "couldn't fit parameters to signature of terra function")
            ctxt:error(ast, retval)
            return errorT
        end
        -- Kinda terrible hack due to flux in Terra inteface here
        if rettype:isstruct() and terralib.sizeof(rettype) == 0 then
            -- case of no return value, no type is needed
            return
        end
        if not T.terraToEbbType(rettype) then
            ctxt:error(ast, "unable to use return type '"..tostring(rettype)..
                            "' of terra function in Ebb")
            return errorT
        end
        return T.terraToEbbType(rettype)
    end
end

function B.terra_to_func(terrafn)
    local newfunc = Builtin.new()
    newfunc.is_a_terra_func = true
    newfunc.terrafn = terrafn
    newfunc.check = TerraCheck(terrafn)
    return newfunc
end
