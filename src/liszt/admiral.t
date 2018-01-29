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

local A = {}
package.loaded["admiral"] = A

import 'ebb'
import 'regent'

local AST  = require 'ebb.src.ast'
local C    = require 'ebb.src.c'
local F    = require 'ebb.src.functions'
local L    = require 'ebblib'
local M    = require 'ebb.src.main'
local P    = require 'ebb.src.phase'
local PRE  = require 'ebb.src.prelude'
local R    = require 'ebb.src.relations'
local RG   = regentlib
local S    = require 'ebb.src.semant'
local T    = require 'ebb.src.types'
local UTIL = require 'ebb.src.util'

local all     = UTIL.all
local isEmpty = UTIL.isEmpty
local keys    = UTIL.keys
local newlist = terralib.newlist

-------------------------------------------------------------------------------
-- Parse config options
-------------------------------------------------------------------------------

local DUMP_REGENT = os.getenv('DUMP_REGENT') == '1'

local LIBS = newlist({'-ljsonparser', '-lm'})

local OBJNAME = os.getenv('OBJNAME') or 'a.out'

local HDF_LIBNAME = os.getenv('HDF_LIBNAME') or 'hdf5'

local HDF_HEADER = os.getenv('HDF_HEADER') or 'hdf5.h'

local USE_HDF = not (os.getenv('USE_HDF') == '0')
if USE_HDF then
  LIBS:insert('-l'..HDF_LIBNAME)
end

-------------------------------------------------------------------------------
-- Helper functions
-------------------------------------------------------------------------------

local function quit(obj)
  print('====================')
  print('Unsupported class:')
  if not obj then
    print('Nil')
  elseif rawget(obj, 'kind') then
    print(obj.kind)
  else
    print('(no kind information)')
  end
  print('--------------------')
  while obj do
    for k,_ in pairs(obj) do
      print(k)
    end
    print('--------------------')
    obj = getmetatable(obj)
  end
  assert(false)
end

-- terralib.type -> number
local function minValue(T)
  return
    (T == int)    and -2147483648 or
    (T == int8)   and -128 or
    (T == int16)  and -32768 or
    (T == int32)  and -2147483648 or
    (T == int64)  and -9223372036854775808 or
    (T == uint)   and 0 or
    (T == uint8)  and 0 or
    (T == uint16) and 0 or
    (T == uint32) and 0 or
    (T == uint64) and 0 or
    (T == bool)   and 0 or
    (T == float)  and -math.huge or
    (T == double) and -math.huge or
    assert(false)
end

-- terralib.type -> number
local function maxValue(T)
  return
    (T == int)    and 2147483647 or
    (T == int8)   and 127 or
    (T == int16)  and 32767 or
    (T == int32)  and 2147483647 or
    (T == int64)  and 9223372036854775807 or
    (T == uint)   and 4294967295 or
    (T == uint8)  and 255 or
    (T == uint16) and 65535 or
    (T == uint32) and 4294967295 or
    (T == uint64) and 18446744073709551615 or
    (T == bool)   and 1 or
    (T == float)  and math.huge or
    (T == double) and math.huge or
    assert(false)
end

-- string, terralib.type -> number
local function opIdentity(op, T)
  return
    (op == '+')   and 0 or
    (op == '-')   and 0 or
    (op == '*')   and 1 or
    (op == '/')   and 1 or
    (op == 'max') and minValue(T) or
    (op == 'min') and maxValue(T) or
    assert(false)
end

-- string -> string
local function opName(op)
  return
    (op == '+')   and 'add' or
    (op == '-')   and 'sub' or
    (op == '*')   and 'mul' or
    (op == '/')   and 'div' or
    (op == 'max') and 'max' or
    (op == 'min') and 'min' or
    assert(false)
end

-- (* -> *) -> ()
local function prettyPrintFun(fun)
  local raw = fun:prettystring(false)
  local startIdx
  for line in raw:gmatch('[^\n]+') do
    if not startIdx then
      local a,b = line:find('%S+')
      startIdx = b + 2
    end
    if not line:find('debuginfo') then
      print(line:sub(startIdx))
    end
  end
end

-- RG.task -> ()
local function prettyPrintTask(tsk)
  tsk:printpretty()
end

-- {string,terralib.type} | {field:string,type:terralib.type} ->
--   string, terralib.type
local function parseStructEntry(entry)
  if terralib.israwlist(entry) and #entry == 2 then
    return entry[1], entry[2]
  elseif entry.field and entry.type then
    return entry.field, entry.type
  else assert(false) end
end

-- terralib.struct -> ()
local function prettyPrintStruct(s)
  print('struct '..s.name..' {')
  for _,e in ipairs(s.entries) do
    local name, type = parseStructEntry(e)
    print('  '..name..' : '..tostring(type)..';')
  end
  print('}')
end

-- string -> string
local function idSanitize(s)
  return s:gsub('[^%w]', '_')
end

local NAME_CACHE = {} -- map(string, (* -> *) | RG.task)

-- RG.task, string -> ()
function A.registerTask(tsk, name)
  name = idSanitize(name)
  while NAME_CACHE[name] do
    name = name..'_'
  end
  NAME_CACHE[name] = tsk
  tsk:set_name(name)
  local parallel_task = tsk:get_parallel_task()
  if parallel_task then
    parallel_task:set_name(name)
  end
  tsk:get_primary_variant():get_ast().name[1] = name -- XXX: Dangerous
  if DUMP_REGENT then
    prettyPrintTask(tsk)
  end
end

-- (* -> *), string -> ()
function A.registerFun(fun, name)
  name = idSanitize(name)
  while NAME_CACHE[name] do
    name = name..'_'
  end
  NAME_CACHE[name] = fun
  fun:setname(name)
  if DUMP_REGENT then
    prettyPrintFun(fun)
  end
end

-- terralib.struct -> ()
function A.registerStruct(s)
  if DUMP_REGENT then
    prettyPrintStruct(s)
  end
end

-------------------------------------------------------------------------------
-- Vector type support
-------------------------------------------------------------------------------

-- terralib.type, int -> terralib.type
local Vector = terralib.memoize(function(T, N)
  return T[N]
end)

-- T:terralib.type, N:int -> (Vector(T,N), Vector(T,N) -> T)
local emitDotProduct = terralib.memoize(function(T, N)
  local a = symbol(Vector(T,N), 'a')
  local b = symbol(Vector(T,N), 'b')
  local elems = newlist()
  local expr = `(a[0] * b[0])
  for i=1,N-1 do
    expr = `expr + (a[ [i] ] * b[ [i] ])
  end
  local terra dot([a], [b]) : T
    return [expr]
  end
  A.registerFun(dot, 'dot_'..tostring(T)..'_'..tostring(N))
  return dot
end)

-- string, T:terralib.type, N:int -> (Vector(T,N), T -> Vector(T,N))
local emitVectorVectorOp = terralib.memoize(function(op, T, N)
  local a = symbol(Vector(T,N), 'a')
  local b = symbol(Vector(T,N), 'b')
  local elems = newlist()
  for i=0,N-1 do
    elems:insert((op == '+') and (`a[ [i] ] + b[ [i] ]) or
                 (op == '-') and (`a[ [i] ] - b[ [i] ]) or
                 (op == '*') and (`a[ [i] ] * b[ [i] ]) or
                 (op == '/') and (`a[ [i] ] / b[ [i] ]) or
                 (op == '%') and (`a[ [i] ] % b[ [i] ]) or
                 assert(false))
  end
  local terra vvop([a], [b]) : Vector(T,N)
    return array([elems])
  end
  A.registerFun(vvop, 'vv_'..opName(op)..'_'..tostring(T)..'_'..tostring(N))
  return vvop
end)

-- string, T:terralib.type, N:int -> (Vector(T,N), T -> Vector(T,N))
local emitVectorScalarOp = terralib.memoize(function(op, T, N)
  local a = symbol(Vector(T,N), 'a')
  local b = symbol(T, 'b')
  local elems = newlist()
  for i=0,N-1 do
    elems:insert((op == '+') and (`a[ [i] ] + b) or
                 (op == '-') and (`a[ [i] ] - b) or
                 (op == '*') and (`a[ [i] ] * b) or
                 (op == '/') and (`a[ [i] ] / b) or
                 (op == '%') and (`a[ [i] ] % b) or
                 assert(false))
  end
  local terra vsop([a], [b]) : Vector(T,N)
    return array([elems])
  end
  A.registerFun(vsop, 'vs_'..opName(op)..'_'..tostring(T)..'_'..tostring(N))
  return vsop
end)

-- terralib.type, RG.rexpr* -> RG.rexpr
local function emitVectorCtor(T, elems)
  return rexpr array([elems]) end
end

-------------------------------------------------------------------------------
-- Translation helper functions
-------------------------------------------------------------------------------

-- () -> bool
function F.Function:isKernel()
  return (#self._decl_ast.params == 1 and
          self._decl_ast.ptypes[1]:iskey() and
          not self._decl_ast.exp)
end

-- map(B.Builtin, (double -> double))
local UNARY_ARITH_FUNS = {
  [L.acos]  = RG.acos(double),
  [L.asin]  = RG.asin(double),
  [L.atan]  = RG.atan(double),
  [L.cbrt]  = RG.cbrt(double),
  [L.ceil]  = RG.ceil(double),
  [L.cos]   = RG.cos(double),
  [L.fabs]  = RG.fabs(double),
  [L.floor] = RG.floor(double),
  [L.log]   = RG.log(double),
  [L.sin]   = RG.sin(double),
  [L.sqrt]  = RG.sqrt(double),
  [L.tan]   = RG.tan(double),
}

-- map(B.Builtin, (double, double -> double))
local BINARY_ARITH_FUNS = {
  [L.pow]   = RG.pow(double),
  [L.fmod]  = RG.fmod(double),
}

-- T.Type -> terralib.type
local function toRType(typ)
  if typ:isprimitive() then
    return typ:terratype()
  elseif typ:isvector() then
    return Vector(toRType(typ.type), typ.N)
  elseif typ:ismatrix() then
    -- TODO: Not supporting matrix types, or operations
    assert(false)
  elseif typ:iskey() then
    return typ.relation:indexType()
  else assert(false) end
end

-- M.ExprConst -> T.Type
local function inferType(lit)
  if type(lit) == 'boolean' then
    return L.bool
  elseif type(lit) == 'number' then
    if lit == math.floor(lit) then
      return L.int
    else
      return L.double
    end
  elseif type(lit) == 'table' then
    assert(terralib.israwlist(lit))
    assert(#lit > 0)
    return L.vector(inferType(lit[1]), #lit)
  else assert(false) end
end

-- M.ExprConst, T.Type? -> RG.rexpr
local function toRConst(lit, typ)
  typ = typ or inferType(lit)
  if typ:iskey() then
    return typ.relation:translateIndex(lit)
  elseif type(lit) == 'boolean' then
    assert(typ == L.bool)
    return rexpr lit end
  elseif type(lit) == 'number' then
    return rexpr [toRType(typ)](lit) end
  elseif type(lit) == 'table' then
    -- TODO: Not supporting matrix types, or operations
    assert(typ:isvector())
    assert(terralib.israwlist(lit))
    assert(#lit == typ.N)
    return emitVectorCtor(
      toRType(typ.type),
      newlist(lit):map(function(e) return toRConst(e, typ.type) end))
  else assert(false) end
end

-- string, T.Type, RG.rexpr, RG.rexpr -> RG.rquote
local function emitReduce(op, typ, lval, exp)
  assert(op ~= '/' or typ == L.float or typ == L.double)
  if typ:isvector() then
    local tmp = RG.newsymbol(toRType(typ), 'tmp')
    local v = RG.newsymbol(toRType(typ), 'v')
    local stmts = newlist()
    stmts:insert(rquote var [tmp] = exp end)
    stmts:insert(rquote var [v] = lval end)
    for i=0,typ.N-1 do
      stmts:insert(
        (op == '+')   and rquote v[ [i] ] +=   tmp[ [i] ]     end or
        (op == '-')   and rquote v[ [i] ] +=   -tmp[ [i] ]    end or
        (op == '*')   and rquote v[ [i] ] *=   tmp[ [i] ]     end or
        (op == '/')   and rquote v[ [i] ] *=   1.0/tmp[ [i] ] end or
        (op == 'max') and rquote v[ [i] ] max= tmp[ [i] ]     end or
        (op == 'min') and rquote v[ [i] ] min= tmp[ [i] ]     end or
        assert(false))
    end
    stmts:insert(rquote lval = v end)
    return rquote [stmts] end
  end
  return
    (op == '+')   and rquote lval +=   exp     end or
    (op == '-')   and rquote lval +=   -exp    end or
    (op == '*')   and rquote lval *=   exp     end or
    (op == '/')   and rquote lval *=   1.0/exp end or
    (op == 'max') and rquote lval max= exp     end or
    (op == 'min') and rquote lval min= exp     end or
    assert(false)
end

-------------------------------------------------------------------------------
-- Global symbols
-------------------------------------------------------------------------------

local NX = RG.newsymbol(int, 'NX')
local NY = RG.newsymbol(int, 'NY')
local NZ = RG.newsymbol(int, 'NZ')

local NUM_PRIM_PARTS = rexpr NX * NY * NZ end

-- () -> RG.symbol
A.primColors = terralib.memoize(function()
  return RG.newsymbol(nil, 'primColors')
end)

-------------------------------------------------------------------------------
-- Relation-to-region translation
-------------------------------------------------------------------------------

-- () -> RG.symbol
R.Relation.regionSymbol = terralib.memoize(function(self)
 return RG.newsymbol(nil, self:Name())
end)

-- () -> RG.symbol
R.Relation.primPartSymbol = terralib.memoize(function(self)
  return RG.newsymbol(nil, self:Name()..'_primPart')
end)

-- () -> RG.symbol
R.Relation.copyRegionSymbol = terralib.memoize(function(self)
 return RG.newsymbol(nil, self:Name()..'_copy')
end)

-- () -> RG.symbol
R.Relation.copyPrimPartSymbol = terralib.memoize(function(self)
  return RG.newsymbol(nil, self:Name()..'_copy_primPart')
end)

-- () -> RG.symbol*
R.Relation.queues = terralib.memoize(function(self)
  assert(self:isCoupled())
  local queues = newlist()
  for i = 1,#self:XferStencil() do
    queues:insert(RG.newsymbol(nil, self:Name()..'_queue_'..(i-1)))
  end
  return queues
end)

-- () -> RG.symbol*
R.Relation.qSrcParts = terralib.memoize(function(self)
  assert(self:isCoupled())
  local qSrcParts = newlist()
  for i = 1,#self:XferStencil() do
    qSrcParts:insert(RG.newsymbol(nil, self:Name()..'_qSrcPart_'..(i-1)))
  end
  return qSrcParts
end)

-- () -> RG.symbol*
R.Relation.qDstParts = terralib.memoize(function(self)
  assert(self:isCoupled())
  local qDstParts = newlist()
  for i = 1,#self:XferStencil() do
    qDstParts:insert(RG.newsymbol(nil, self:Name()..'_qDstPart_'..(i-1)))
  end
  return qDstParts
end)

-- () -> RG.index_type
function R.Relation:indexType()
  if self:isPlain() or self:isCoupled() then
    return int1d
  elseif self:isGrid() then
    return int3d
  else assert(false) end
end

-- () -> RG.ispace_type
function R.Relation:indexSpaceType()
  return ispace(self:indexType())
end

-- () -> terralib.struct
R.Relation.fieldSpace = terralib.memoize(function(self)
  local fs = terralib.types.newstruct(self:Name() .. '_columns')
  for _,fld in ipairs(self:Fields()) do
    fs.entries:insert({field=fld:Name(), type=toRType(fld:Type())})
  end
  if self:isCoupled() then
    fs.entries:insert({field='__valid', type=bool})
  end
  A.registerStruct(fs)
  return fs
end)

-- () -> RG.region_type
function R.Relation:regionType()
  -- Region types in signatures must be distinct, so we're not caching here.
  return region(self:indexSpaceType(), self:fieldSpace())
end

-- () -> terralib.array
R.Relation.queueFieldSpace = terralib.memoize(function(self)
  assert(self:isCoupled())
  return int8[terralib.sizeof(self:fieldSpace())]
end)

-- () -> RG.region_type
function R.Relation:queueRegionType()
  assert(self:isCoupled())
  -- Region types in signatures must be distinct, so we're not caching here.
  return region(self:indexSpaceType(), self:queueFieldSpace())
end

-- () -> int
R.Relation.validFieldOffset = terralib.memoize(function(self)
  assert(self:isCoupled())
  -- Ask Terra the offset to the __valid field
  local terra getOffset()
    var x : self:fieldSpace()
    return [int64]([&int8](&(x.__valid)) - [&int8](&x))
  end
  A.registerFun(getOffset, self:Name()..'_getOffset')
  return getOffset()
end)

-- M.ExprConst -> RG.rexpr
function R.Relation:translateIndex(lit)
  if self:isPlain() or self:isCoupled() then
    return rexpr [int1d]([toRConst(lit, L.int)]) end
  elseif self:isGrid() then
    assert(terralib.israwlist(lit))
    assert(#lit == 3)
    return rexpr [int3d]({ x = [toRConst(lit[1], L.int)],
                           y = [toRConst(lit[2], L.int)],
                           z = [toRConst(lit[3], L.int)]}) end
  else assert(false) end
end

-- () -> RG.rexpr
function R.Relation:primPartSize()
  assert(self:isCoupled())
  return rexpr C.ceil([self:Size():varSymbol()] / NUM_PRIM_PARTS
                      * [self:MaxSkew():varSymbol()]) end
end

-- () -> RG.rexpr
function R.Relation:emitISpaceInit()
  if self:isPlain() then
    return rexpr ispace(int1d, [self:Size():varSymbol()]) end
  elseif self:isGrid() then
    local xNum = self:xNum():varSymbol()
    local yNum = self:yNum():varSymbol()
    local zNum = self:zNum():varSymbol()
    local xBnum = self:xBnum():varSymbol()
    local yBnum = self:yBnum():varSymbol()
    local zBnum = self:zBnum():varSymbol()
    return rexpr
      ispace(int3d, { x = xNum + 2 * xBnum,
                      y = yNum + 2 * yBnum,
                      z = zNum + 2 * zBnum })
    end
  elseif self:isCoupled() then
    return rexpr ispace(int1d, [self:primPartSize()] * NUM_PRIM_PARTS) end
  else assert(false) end
end

-- () -> RG.rquote
function R.Relation:emitRegionInit()
  local fs = self:fieldSpace()
  return rquote
    var is = [self:emitISpaceInit()]
    var [self:regionSymbol()] = region(is, fs)
    -- TODO:
    -- single copy-region for every dump/load variant
    -- declaring all fields on the copy-region, not just the ones dumped/loaded
    var [self:copyRegionSymbol()] = region(is, fs)
  end
end

-- () -> RG.rquote
function R.Relation:emitValidInit()
  assert(self:isCoupled())
  local initValidField
  __demand(__parallel) task initValidField(r : self:regionType())
  where writes(r.__valid) do
    for e in r do
      e.__valid = false
    end
  end
  A.registerTask(initValidField, self:Name()..'_initValidField')
  return rquote initValidField([self:regionSymbol()]) end
end

-- () -> RG.rquote
function R.Relation:emitPrimPartInit()
  local r = self:regionSymbol()
  local p_r = self:primPartSymbol()
  local s = self:copyRegionSymbol()
  local p_s = self:copyPrimPartSymbol()
  if self:isPlain() then
    return rquote
      RG.assert([self:Size():varSymbol()] % NUM_PRIM_PARTS == 0,
                'Uneven partitioning')
      var [p_r] = partition(equal, r, [A.primColors()])
      var [p_s] = partition(equal, s, [A.primColors()])
    end
  elseif self:isGrid() then
    local xNum = self:xNum():varSymbol()
    local yNum = self:yNum():varSymbol()
    local zNum = self:zNum():varSymbol()
    local xBnum = self:xBnum():varSymbol()
    local yBnum = self:yBnum():varSymbol()
    local zBnum = self:zBnum():varSymbol()
    return rquote
      RG.assert(xNum % NX == 0, 'Uneven partitioning')
      RG.assert(yNum % NY == 0, 'Uneven partitioning')
      RG.assert(zNum % NZ == 0, 'Uneven partitioning')
      var coloring = RG.c.legion_domain_point_coloring_create()
      -- Partition interior equally, then add boundaries
      for c in [A.primColors()] do
        var rect = rect3d{
          lo = int3d{ x = xBnum + (xNum / NX) * c.x,
                      y = yBnum + (yNum / NY) * c.y,
                      z = zBnum + (zNum / NZ) * c.z },
          hi = int3d{ x = xBnum + (xNum / NX) * (c.x+1) - 1,
                      y = yBnum + (yNum / NY) * (c.y+1) - 1,
                      z = zBnum + (zNum / NZ) * (c.z+1) - 1 }
        }
        if c.x == 0    then rect.lo.x -= xBnum end
        if c.x == NX-1 then rect.hi.x += xBnum end
        if c.y == 0    then rect.lo.y -= yBnum end
        if c.y == NY-1 then rect.hi.y += yBnum end
        if c.z == 0    then rect.lo.z -= zBnum end
        if c.z == NZ-1 then rect.hi.z += zBnum end
        RG.c.legion_domain_point_coloring_color_domain(coloring, c, rect)
      end
      var [p_r] = partition(disjoint, r, coloring, [A.primColors()])
      var [p_s] = partition(disjoint, s, coloring, [A.primColors()])
      RG.c.legion_domain_point_coloring_destroy(coloring)
    end
  elseif self:isCoupled() then
    local primPartSize = self:primPartSize()
    return rquote
      RG.assert([self:Size():varSymbol()] % NUM_PRIM_PARTS == 0,
                'Uneven partitioning')
      var coloring = RG.c.legion_domain_point_coloring_create()
      for z = 0, NZ do
        for y = 0, NY do
          for x = 0, NX do
            var rBase : int64
            for rStart in r do
              rBase = rStart + (z*NX*NY + y*NX + x) * primPartSize
              break
            end
            RG.c.legion_domain_point_coloring_color_domain(
              coloring, int3d{x,y,z},
              [rect1d]{ rBase, rBase + primPartSize - 1 })
          end
        end
      end
      var [p_r] = partition(disjoint, r, coloring, [A.primColors()])
      var [p_s] = partition(disjoint, s, coloring, [A.primColors()])
      RG.c.legion_domain_point_coloring_destroy(coloring)
    end
  else assert(false) end
end

-- () -> RG.task
R.Relation.emitElemColor = terralib.memoize(function(self)
  local NX_ = RG.newsymbol(int, 'NX_')
  local NY_ = RG.newsymbol(int, 'NY_')
  local NZ_ = RG.newsymbol(int, 'NZ_')
  local elemColor
  if self:isGrid() then
    __demand(__inline) task elemColor(idx : int3d,
                                      xNum : int, yNum : int, zNum : int,
                                      xBnum : int, yBnum : int, zBnum : int,
                                      [NX_], [NY_], [NZ_])
      idx.x = min( max( idx.x, xBnum ), xNum + xBnum - 1 )
      idx.y = min( max( idx.y, yBnum ), yNum + yBnum - 1 )
      idx.z = min( max( idx.z, zBnum ), zNum + zBnum - 1 )
      return int3d{ (idx.x - xBnum) / (xNum / NX_),
                    (idx.y - yBnum) / (yNum / NY_),
                    (idx.z - zBnum) / (zNum / NZ_) }
    end
  else
    -- TODO: Not covered: plain and coupled relations. These follow a
    -- different partitioning scheme, and may also require the base region to
    -- calculate offsets.
    assert(false)
  end
  A.registerTask(elemColor, self:Name()..'_elemColor')
  return elemColor
end)

-- int -> RG.rquote
function R.Relation:emitQueueInit(i)
  assert(self:isCoupled())
  local q = self:queues()[i]
  local fspaceExpr = self:queueFieldSpace()
  local maxXferNum = self:MaxXferNum():varSymbol()
  return rquote
    var [q] = region(ispace(int1d, maxXferNum * NUM_PRIM_PARTS), fspaceExpr)
  end
end

-- int -> RG.rquote
function R.Relation:emitQueuePartInit(i)
  assert(self:isCoupled())
  local q = self:queues()[i]
  local qSrcPart = self:qSrcParts()[i]
  local qDstPart = self:qDstParts()[i]
  local colors = A.primColors()
  local stencil = self:XferStencil()[i]
  local maxXferNum = self:MaxXferNum():varSymbol()
  return rquote
    var srcColoring = RG.c.legion_domain_point_coloring_create()
    for z = 0, NZ do
      for y = 0, NY do
        for x = 0, NX do
          var qBase : int64
          for qStart in q do
            qBase = qStart + (z*NX*NY + y*NX + x) * maxXferNum
            break
          end
          RG.c.legion_domain_point_coloring_color_domain(
            srcColoring, int3d{x,y,z},
            [rect1d]{ qBase, qBase + maxXferNum - 1 })
        end
      end
    end
    var [qSrcPart] = partition(disjoint, q, srcColoring, colors)
    RG.c.legion_domain_point_coloring_destroy(srcColoring)
    var dstColoring = RG.c.legion_domain_point_coloring_create()
    var colorOff = int3d{ [stencil[1]], [stencil[2]], [stencil[3]] }
    for c in colors do
      var srcBase : int64
      for qptr in qSrcPart[(c - colorOff + {NX,NY,NZ}) % {NX,NY,NZ}] do
        srcBase = [int1d](qptr)
        break
      end
      RG.c.legion_domain_point_coloring_color_domain(
        dstColoring, c, [rect1d]{ srcBase, srcBase + maxXferNum - 1 })
    end
    var [qDstPart] = partition(aliased, q, dstColoring, colors)
    RG.c.legion_domain_point_coloring_destroy(dstColoring)
  end
end

-- () -> RG.task
R.Relation.emitPushAll = terralib.memoize(function(self)
  assert(self:isCoupled())
  local couplingFld = self:CouplingField()
  local NX_ = RG.newsymbol(int, 'NX_')
  local NY_ = RG.newsymbol(int, 'NY_')
  local NZ_ = RG.newsymbol(int, 'NZ_')
  -- utilized sub-tasks
  local rngElemColor = couplingFld:Type().relation:emitElemColor()
  -- code elements to fill in
  local qs = newlist() -- RG.symbol*
  local privileges = newlist() -- RG.privilege*
  local queueInits = newlist() -- RG.rquote*
  local moveChecks = newlist() -- RG.rquote*
  -- visible symbols
  local r = RG.newsymbol(self:regionType(), 'r')
  local rPtr = RG.newsymbol(nil, 'rPtr')
  local partColor = RG.newsymbol(int3d, 'partColor')
  local elemColor = RG.newsymbol(int3d, 'elemColor')
  -- fill in parameters & movement cases
  local terra pushElement(dst : &opaque,
                          idx : int,
                          src : self:regionType().fspace_type)
    var ptr : &int8 = [&int8](dst) + idx * [self:queueFieldSpace().N]
    C.memcpy(ptr, &src, [self:queueFieldSpace().N])
  end
  A.registerFun(pushElement, self:Name()..'_pushElement')
  local terra getBasePointer(pr : RG.c.legion_physical_region_t,
                             fid : RG.c.legion_field_id_t,
                             runtime : RG.c.legion_runtime_t)
    var acc = RG.c.legion_physical_region_get_field_accessor_array_1d(pr, fid)
    var lr = RG.c.legion_physical_region_get_logical_region(pr)
    var domain = RG.c.legion_index_space_get_domain(runtime, lr.index_space)
    var rect = RG.c.legion_domain_get_rect_1d(domain)
    var subrect : RG.c.legion_rect_1d_t
    var offsets : RG.c.legion_byte_offset_t[1]
    var p = RG.c.legion_accessor_array_1d_raw_rect_ptr(
        acc, rect, &subrect, &(offsets[0]))
    RG.c.legion_accessor_array_1d_destroy(acc)
    return p
  end
  A.registerFun(getBasePointer, self:Name()..'_getBasePointer')
  for i,stencil in ipairs(self:XferStencil()) do
    local q = RG.newsymbol(self:queueRegionType(), 'q'..(i-1))
    local qBasePtr = RG.newsymbol(&opaque, 'qBasePtr'..(i-1))
    qs:insert(q)
    privileges:insert(RG.privilege(RG.reads, q))
    privileges:insert(RG.privilege(RG.writes, q))
    queueInits:insert(rquote
      for qPtr in [q] do
        [q][qPtr][ [self:validFieldOffset()] ] = [int8](false)
      end
      var [qBasePtr] =
        getBasePointer(__physical([q])[0], __fields([q])[0], __runtime())
    end)
    moveChecks:insert(rquote
      do
        var colorOff = int3d{ [stencil[1]], [stencil[2]], [stencil[3]] }
        if rPtr.__valid and
           elemColor == (partColor + colorOff + {NX_,NY_,NZ_}) % {NX_,NY_,NZ_} then
          var idx = 0
          for qPtr in [q] do
            if not [bool]([q][qPtr][ [self:validFieldOffset()] ]) then
              pushElement(qBasePtr, idx, r[rPtr])
              rPtr.__valid = false
              RG.assert([bool]([q][qPtr][ [self:validFieldOffset()] ]),
                        'Element did not get copied properly')
              break
            end
            idx += 1
          end
          RG.assert(not rPtr.__valid, 'Transfer queue ran out of space')
        end
      end
    end)
  end
  -- synthesize task
  local task pushAll([partColor], [r], [qs],
                     rngXNum : int, rngYNum : int, rngZNum : int,
                     rngXbnum : int, rngYbnum : int, rngZbnum : int,
                     [NX_], [NY_], [NZ_])
  where reads(r), writes(r.__valid), [privileges] do
    [queueInits]
    for [rPtr] in r do
      if rPtr.__valid then
        var [elemColor] =
          rngElemColor(rPtr.[couplingFld:Name()],
                       rngXNum, rngYNum, rngZNum,
                       rngXbnum, rngYbnum, rngZbnum,
                       NX_, NY_, NZ_)
        if elemColor ~= partColor then
          [moveChecks]
          RG.assert(not rPtr.__valid, 'Element moved past predicted stencil')
        end
      end
    end
  end
  A.registerTask(pushAll, self:Name()..'_pushAll')
  return pushAll
end)

-- () -> RG.task
R.Relation.emitPullAll = terralib.memoize(function(self)
  assert(self:isCoupled())
  local queues = newlist() -- RG.symbol*
  local privileges = newlist() -- RG.privilege*
  for i = 1,#self:XferStencil() do
    local q = RG.newsymbol(self:queueRegionType(), 'q'..(i-1))
    queues:insert(q)
    privileges:insert(RG.privilege(RG.reads, q))
  end
  local terra pullElement(src : &self:queueRegionType().fspace_type.type)
    var dst : self:regionType().fspace_type
    C.memcpy(&dst, src, [self:queueRegionType().fspace_type.N])
    return dst
  end
  A.registerFun(pullElement, self:Name()..'_pullElement')
  local task pullAll(color : int3d, r : self:regionType(), [queues])
  where reads writes(r), [privileges] do
    [queues:map(function(q) return rquote
      for qPtr in q do
        -- TODO: Check that the element is coming to the appropriate partition.
        if [bool](q[qPtr][ [self:validFieldOffset()] ]) then
          var copied = false
          for rPtr in r do
            if not rPtr.__valid then
              r[rPtr] = pullElement(q[qPtr])
              copied = true
              RG.assert(r[rPtr].__valid,
                'Pulled particle was not copied correctly')
              break
            end
          end
          RG.assert(copied, 'Ran out of space on sub-partition')
        end
      end
    end end)]
  end
  A.registerTask(pullAll, self:Name()..'_pullAll')
  return pullAll
end)

-- () -> RG.rquote
function R.Relation:emitPrimPartUpdate()
  assert(self:isCoupled())
  local pushAll = self:emitPushAll()
  local pullAll = self:emitPullAll()
  local rngRel = self:CouplingField():Type().relation
  local rngXNum = rngRel:xNum():varSymbol()
  local rngYNum = rngRel:yNum():varSymbol()
  local rngZNum = rngRel:zNum():varSymbol()
  local rngXbnum = rngRel:xBnum():varSymbol()
  local rngYbnum = rngRel:yBnum():varSymbol()
  local rngZbnum = rngRel:zBnum():varSymbol()
  return rquote
    for c in [A.primColors()] do
      pushAll(c, [self:primPartSymbol()][c],
              [self:qSrcParts():map(function(qSrcPart)
                return rexpr qSrcPart[c] end
              end)],
              rngXNum, rngYNum, rngZNum,
              rngXbnum, rngYbnum, rngZbnum,
              NX, NY, NZ)
    end
    for c in [A.primColors()] do
      pullAll(c, [self:primPartSymbol()][c],
              [self:qDstParts():map(function(qDstPart)
                return rexpr qDstPart[c] end
              end)])
    end
    -- TODO: Check that all out-queues have been emptied out.
  end
end

-- () -> RG.rquote
function R.Relation:emitPrimPartCheck()
  assert(self:isCoupled())
  local couplingFld = self:CouplingField()
  local domPart = self:primPartSymbol()
  local rngPart = couplingFld:Type().relation:primPartSymbol()
  return rquote
    for c in [A.primColors()] do
      for x in [domPart][c] do
        if x.__valid then
          -- TODO: This membership check is inefficient; we do the original
          -- partitioning manually, so we know the exact limits of the
          -- sub-partitions.
          var found = false
          for y in [rngPart][c] do
            if x.[couplingFld:Name()] == y then
              found = true
              break
            end
          end
          RG.assert(found, 'Automatic partitioning invariant violated')
        end
      end
    end
  end
end

-------------------------------------------------------------------------------
-- Function translation
-------------------------------------------------------------------------------

-- FunInfo = {
--   name           : string
--   domainRel      : R.Relation?
--   field_use      : map(R.Field, P.PhaseType)
--   global_use     : map(PRE.Global, P.PhaseType)
--   inserts        : map(R.Relation, AST.InsertStatement)
--   deletes        : map(R.Relation, AST.DeleteStatement)
--   manualParal    : bool
-- }

local FunContext = {}
FunContext.__index = FunContext

-- FunInfo, AST.Symbol*, T.Type* -> FunContext
function FunContext.New(info, argNames, argTypes)
  local self = setmetatable({
    -- Symbol mappings
    localMap        = {},        -- map(AST.Symbol, RG.rexpr)
    globalMap       = {},        -- map(PRE.Global, RG.symbol)
    relMap          = {},        -- map(R.Relation, RG.symbol)
    -- Signature information
    domainRel       = nil,       -- R.Relation?
    args            = newlist(), -- RG.Symbol*
    accessedRels    = newlist(), -- R.Relation*
    readGlobals     = newlist(), -- PRE.Global*
    -- Field use information
    privileges      = newlist(), -- RG.privilege*
    -- Global reduction information
    reducedGlobal   = nil,       -- PRE.Global?
    globalReduceAcc = nil,       -- RG.symbol?
    globalReduceOp  = nil,       -- string?
  }, FunContext)
  -- Process mapped relation
  if info.domainRel then
    local rel = info.domainRel
    self.domainRel = rel
    -- Mapped relation always set as first accessed relation.
    self.accessedRels:insert(rel)
    self.relMap[rel] = RG.newsymbol(rel:regionType(), rel:Name())
  end
  -- Process arguments
  for i,lsym in ipairs(argNames) do
    local T
    if i > 1 or not info.domainRel then
      -- If this is a kernel, leave the first argument untyped (it will be used
      -- as the loop variable, and not included in the signature).
      T = toRType(argTypes[i])
    end
    local rsym = RG.newsymbol(T, tostring(lsym))
    self.args:insert(rsym)
    self.localMap[lsym] = rsym
  end
  -- Process field access modes
  local fields = keys(info.field_use)
  fields:sort(function(f1,f2) return f1:FullName() < f2:FullName() end)
  for _,fld in ipairs(fields) do
    local pt = info.field_use[fld]
    local rel = fld:Relation()
    local rg = self.relMap[rel]
    if not rg then
      rg = RG.newsymbol(rel:regionType(), rel:Name())
      self.accessedRels:insert(rel)
      self.relMap[rel] = rg
    end
    -- Assuming phase checking has checked for errors
    if pt.read or pt.write then
      self.privileges:insert(RG.privilege(RG.reads, rg, fld:Name()))
    end
    if pt.write then
      self.privileges:insert(RG.privilege(RG.writes, rg, fld:Name()))
    end
    if pt.reduceop then
      self.privileges:insert(RG.privilege(RG.reads, rg, fld:Name()))
      self.privileges:insert(RG.privilege(RG.writes, rg, fld:Name()))
    end
  end
  -- Process inserts and deletes
  for rel,_ in pairs(info.inserts or {}) do
    assert(not self.relMap[rel])
    local rg = RG.newsymbol(rel:regionType(), rel:Name())
    self.accessedRels:insert(rel)
    self.relMap[rel] = rg
    for _,fld in ipairs(rel:Fields()) do
      self.privileges:insert(RG.privilege(RG.reads, rg, fld:Name()))
      self.privileges:insert(RG.privilege(RG.writes, rg, fld:Name()))
    end
    self.privileges:insert(RG.privilege(RG.reads, rg, '__valid'))
    self.privileges:insert(RG.privilege(RG.writes, rg, '__valid'))
  end
  for rel,_ in pairs(info.deletes or {}) do
    local rg = assert(self.relMap[rel])
    self.privileges:insert(RG.privilege(RG.reads, rg, '__valid'))
    self.privileges:insert(RG.privilege(RG.writes, rg, '__valid'))
  end
  -- Privileges for accessing translator-added flags
  if info.domainRel and info.domainRel:isCoupled() then
    local rg = assert(self.relMap[info.domainRel])
    self.privileges:insert(RG.privilege(RG.reads, rg, '__valid'))
  end
  -- Process global access modes
  local globals = keys(info.global_use)
  globals:sort(function(g1,g2) return g1:Name() < g2:Name() end)
  for _,g in ipairs(globals) do
    local pt = info.global_use[g]
    if pt.read and not pt.reduceop then
      assert(not self.globalMap[g])
      self.globalMap[g] = RG.newsymbol(toRType(g:Type()), idSanitize(g:Name()))
      self.readGlobals:insert(g)
    elseif pt.reduceop and not pt.read then
      assert(not self.reducedGlobal)
      self.reducedGlobal = g
      self.globalReduceAcc = RG.newsymbol(toRType(g:Type()), 'acc')
      self.globalReduceOp = pt.reduceop
    else assert(false) end
  end
  return self
end

-- PRE.Global -> RG.symbol
function FunContext:addGlobal(global)
  local sym = self.globalMap[global]
  if not sym then
    self.readGlobals:insert(global)
    sym = RG.newsymbol(toRType(global:Type()), idSanitize(global:Name()))
    self.globalMap[global] = sym
  end
  return sym
end

-- AST.Symbol, T.Type -> RG.symbol
function FunContext:addLocal(lsym, type)
  assert(not self.localMap[lsym])
  local rsym = RG.newsymbol(toRType(type), tostring(lsym))
  self.localMap[lsym] = rsym
  return rsym
end

-- AST.Symbol, RG.rexpr -> ()
function FunContext:addAlias(lsym, expr)
  assert(not self.localMap[lsym])
  self.localMap[lsym] = expr
end

-- () -> RG.symbol*
function FunContext:signature()
  local fullArgs = newlist()
  for i = 1, #self.args do
    if i > 1 or not self.domainRel then
      fullArgs:insert(self.args[i])
    end
  end
  for i = 1, #self.accessedRels do
    local rel = self.accessedRels[i]
    fullArgs:insert(self.relMap[rel])
  end
  for i = 1, #self.readGlobals do
    local g = self.readGlobals[i]
    fullArgs:insert(self.globalMap[g])
  end
  return fullArgs
end

-- () -> bool
function AST.Expression:hasNoStencil()
  error('Abstract Method')
end
function AST.BinaryOp:hasNoStencil()
  return self.lhs:hasNoStencil() and self.rhs:hasNoStencil()
end
function AST.Bool:hasNoStencil()
  return true
end

function AST.Call:hasNoStencil()
  if self.func == L.print or self.func == L.rand or
     self.func == L.id or
     self.func == L.xid or self.func == L.yid or self.func == L.zid or
     UNARY_ARITH_FUNS[self.func] or BINARY_ARITH_FUNS[self.func] or
     self.func == L.fmax or self.func == L.imax or
     self.func == L.fmin or self.func == L.imin or
     self.func == L.dot or self.func == L.times then
    return all(self.params, function(e) return e:hasNoStencil() end)
  end
  if self.func == L.UNSAFE_ROW then
    return self.params[1]:hasNoStencil()
  end
  return false
end
function AST.Cast:hasNoStencil()
  return self.value:hasNoStencil()
end
function AST.FieldAccess:hasNoStencil()
  return self.key:hasNoStencil()
end
function AST.FieldAccessIndex:hasNoStencil()
  return self.base:hasNoStencil() and self.index:hasNoStencil()
end
function AST.Global:hasNoStencil()
  return true
end
function AST.GlobalIndex:hasNoStencil()
  return self.index:hasNoStencil()
end
function AST.LetExpr:hasNoStencil()
  return self.block:hasNoStencil() and self.exp:hasNoStencil()
end
function AST.LuaObject:hasNoStencil()
  quit(self)
end
function AST.MatrixLiteral:hasNoStencil()
  quit(self)
end
function AST.Name:hasNoStencil()
  return true
end
function AST.Number:hasNoStencil()
  return true
end
function AST.Quote:hasNoStencil()
  return self.code:hasNoStencil()
end
function AST.RecordLiteral:hasNoStencil()
  quit(self)
end
function AST.Reduce:hasNoStencil()
  quit(self)
end
function AST.SquareIndex:hasNoStencil()
  return self.base:hasNoStencil() and self.index:hasNoStencil()
end
function AST.String:hasNoStencil()
  quit(self)
end
function AST.TableLookup:hasNoStencil()
  quit(self)
end
function AST.UnaryOp:hasNoStencil()
  return self.exp:hasNoStencil()
end
function AST.VectorLiteral:hasNoStencil()
  return all(self.elems, function(e) return e:hasNoStencil() end)
end
function AST.Where:hasNoStencil()
  quit(self)
end

-- () -> bool
function AST.Statement:hasNoStencil()
  error('Abstract Method')
end
function AST.Assignment:hasNoStencil()
  return self.lvalue:hasNoStencil() and self.exp:hasNoStencil()
end
function AST.Break:hasNoStencil()
  quit(self)
end
function AST.DeclStatement:hasNoStencil()
  return not self.initializer or self.initializer:hasNoStencil()
end
function AST.DeleteStatement:hasNoStencil()
  return self.key:hasNoStencil()
end
function AST.DoStatement:hasNoStencil()
  return self.body:hasNoStencil()
end
function AST.ExprStatement:hasNoStencil()
  return self.exp:hasNoStencil()
end
function AST.FieldWrite:hasNoStencil()
  return self.fieldaccess:hasNoStencil() and self.exp:hasNoStencil()
end
function AST.GenericFor:hasNoStencil()
  quit(self)
end
function AST.GlobalReduce:hasNoStencil()
  return self.exp:hasNoStencil()
end
function AST.IfStatement:hasNoStencil()
  return all(self.if_blocks, function(cb)
               return cb.cond:hasNoStencil() and cb.body:hasNoStencil()
             end)
    and (not self.else_block or self.else_block:hasNoStencil())
end
function AST.InsertStatement:hasNoStencil()
  return all(self.record.exprs, function(e) return e:hasNoStencil() end)
end
function AST.NumericFor:hasNoStencil()
  return self.body:hasNoStencil()
end
function AST.RepeatStatement:hasNoStencil()
  quit(self)
end
function AST.WhileStatement:hasNoStencil()
  quit(self)
end

-- () -> bool
function AST.Block:hasNoStencil()
  return all(self.statements, function(s) return s:hasNoStencil() end)
end

-- FunInfo -> bool
local function mustParallelizeManually(info)
  -- HACK: Need to manually parallelize kernels that insert or delete
  return not (isEmpty(info.inserts) and isEmpty(info.deletes))
         or info.manualParal
end

-- FunInfo -> RG.task, FunContext
function AST.UserFunction:toTask(info)
  -- self.params : AST.Symbol*
  -- self.ptypes : T.Type*
  -- self.body   : AST.Block
  -- self.exp    : AST.Expression?
  local ctxt = FunContext.New(info, self.params, self.ptypes)
  assert(not ctxt.reducedGlobal or not self.exp)
  -- Synthesize body
  local body = newlist() -- RG.quote*
  if ctxt.reducedGlobal then
    local accInitVal =
      opIdentity(ctxt.globalReduceOp, ctxt.globalReduceAcc:gettype())
    body:insert(rquote var [ctxt.globalReduceAcc] = [accInitVal] end)
  end
  local block = self.body:toRQuote(ctxt)
  if info.domainRel then
    local loopVar = ctxt.args[1]
    local rel = ctxt.relMap[info.domainRel]
    if info.domainRel:isCoupled() then
      block = rquote if [rel][loopVar].__valid then [block] end end
    end
    if RG.config['openmp'] then
      block = rquote
        __demand(__openmp) for [loopVar] in [rel] do [block] end
      end
    else
      block = rquote for [loopVar] in [rel] do [block] end end
    end
  end
  body:insert(block)
  if ctxt.reducedGlobal then
    body:insert(rquote return [ctxt.globalReduceAcc] end)
  end
  if self.exp then
    body:insert(rquote return [self.exp:toRExpr(ctxt)] end)
  end
  -- Synthesize task
  local tsk
  if info.domainRel then
    if mustParallelizeManually(info) then
      -- TODO: Only handling the simple case of functions without stencils,
      -- which don't require any changes.
      assert(self.body:hasNoStencil())
      task tsk([ctxt:signature()])
      where [ctxt.privileges] do [body] end
    elseif not ctxt.reducedGlobal and RG.check_cuda_available() then
      __demand(__parallel, __cuda) task tsk([ctxt:signature()])
      where [ctxt.privileges] do [body] end
    else
      __demand(__parallel) task tsk([ctxt:signature()])
      where [ctxt.privileges] do [body] end
    end
  else
    __demand(__inline) task tsk([ctxt:signature()])
    where [ctxt.privileges] do [body] end
  end
  -- Finalize task
  A.registerTask(tsk, info.name)
  return tsk, ctxt
end

local TO_KERNEL_TASK_CACHE = {} -- map(F.Function, {RG.task,FunContext})

-- () -> RG.task, FunContext
function F.Function:toKernelTask()
  if TO_KERNEL_TASK_CACHE[self] then
    return unpack(TO_KERNEL_TASK_CACHE[self])
  end
  assert(self:isKernel())
  local argRel = self._decl_ast.ptypes[1].relation
  local info = self:_get_typechecked(42, argRel, {})
  local typedAST = info.typed_ast
  -- info : {
  --   ...
  --   typed_ast      : AST.UserFunction
  --   field_use      : map(R.Field, P.PhaseType)
  --   global_use     : map(PRE.Global, P.PhaseType)
  --   inserts        : map(R.Relation, AST.InsertStatement)
  --   deletes        : map(R.Relation, AST.DeleteStatement)
  -- }
  info.name = self:Name()
  info.domainRel = argRel
  info.manualParal = self._MANUAL_PARAL
  local tsk, ctxt = typedAST:toTask(info)
  TO_KERNEL_TASK_CACHE[self] = {tsk, ctxt}
  return tsk, ctxt
end

local TO_HELPER_TASK_CACHE = {} -- map(F.Function, {RG.task,FunContext})

-- T.Type* -> RG.task, FunContext
function F.Function:toHelperTask(argTypes)
  -- TODO: Only caching on the function object; we assume the helper functions
  -- have a single specialization.
  if TO_HELPER_TASK_CACHE[self] then
    return unpack(TO_HELPER_TASK_CACHE[self])
  end
  local typedAST = S.check_helper_func(self, argTypes)
  local info = P.phasePass(typedAST)
  -- info : {
  --   ...
  --   field_use      : map(R.Field, P.PhaseType)
  --   global_use     : map(PRE.Global, P.PhaseType)
  --   inserts        : map(R.Relation, AST.InsertStatement)
  --   deletes        : map(R.Relation, AST.DeleteStatement)
  -- }
  info.name = self:Name()
  info.domainRel = nil
  info.manualParal = self._MANUAL_PARAL
  typedAST.ptypes = argTypes
  local tsk, ctxt = typedAST:toTask(info)
  TO_HELPER_TASK_CACHE[self] = {tsk, ctxt}
  return tsk, ctxt
end

-------------------------------------------------------------------------------
-- AST translation
-------------------------------------------------------------------------------

-- AST.Expression, FunContext -> RG.rexpr
local function recoverHelperCall(expr, ctxt)
  -- expr.orig_func   : F.Function
  -- expr.orig_params : AST.Expression*
  local argTypes = newlist(expr.orig_params):map(function(p)
    return p.node_type
  end)
  local hTask, hCtxt = expr.orig_func:toHelperTask(argTypes)
  local actualArgs = newlist()
  for i = 1, #expr.orig_params do
    local p = expr.orig_params[i]
    actualArgs:insert(p:toRExpr(ctxt))
  end
  for i = 1, #hCtxt.accessedRels do
    local rel = hCtxt.accessedRels[i]
    actualArgs:insert(assert(ctxt.relMap[rel]))
  end
  for i = 1, #hCtxt.readGlobals do
    local g = hCtxt.readGlobals[i]
    actualArgs:insert(assert(ctxt.globalMap[g]))
  end
  return rexpr [hTask]([actualArgs]) end
end

-- FunContext -> RG.rexpr
function AST.Expression:toRExpr(ctxt)
  error('Abstract Method')
end
function AST.BinaryOp:toRExpr(ctxt)
  -- self.lhs : AST.Expression
  -- self.rhs : AST.Expression
  -- self.op  : string
  local a = self.lhs:toRExpr(ctxt)
  local b = self.rhs:toRExpr(ctxt)
  local t1 = self.lhs.node_type
  local t2 = self.rhs.node_type
  if t1:isvector() and t2:isvector() then
    local fun = emitVectorVectorOp(self.op, toRType(t1.type), t1.N)
    return rexpr fun(a, b) end
  elseif t1:isvector() and t2:isscalar() then
    local fun = emitVectorScalarOp(self.op, toRType(t2), t1.N)
    return rexpr fun(a, b) end
  elseif t1:isscalar() and t2:isvector() and self.op == '*' then
    local fun = emitVectorScalarOp(self.op, toRType(t1), t2.N)
    return rexpr fun(b, a) end
  end
  assert(t1:isscalar() and t2:isscalar())
  return
    (self.op == '==')  and rexpr a == b  end or
    (self.op == '~=')  and rexpr a ~= b  end or
    (self.op == '<')   and rexpr a < b   end or
    (self.op == '>')   and rexpr a > b   end or
    (self.op == '<=')  and rexpr a <= b  end or
    (self.op == '>=')  and rexpr a >= b  end or
    (self.op == 'or')  and rexpr a or b  end or
    (self.op == 'and') and rexpr a and b end or
    (self.op == '+')   and rexpr a + b   end or
    (self.op == '-')   and rexpr a - b   end or
    (self.op == '*')   and rexpr a * b   end or
    (self.op == '/')   and rexpr a / b   end or
    (self.op == '%')   and rexpr a % b   end or
    assert(false)
end
function AST.Bool:toRExpr(ctxt)
  -- self.value : bool
  return rexpr [self.value] end
end
function AST.Call:toRExpr(ctxt)
  -- self.func   : B.Builtin
  -- self.params : table*
  assert(L.is_builtin(self.func))
  -- Affine expression
  -- self.params[1] : AST.LuaObject
  --   .node_type.value : R.Relation
  -- self.params[2] : AST.MatrixLiteral
  --   .n       : int
  --   .m       : int
  --   .matvals : int[.n][.m]
  -- self.params[3] : AST.Expression
  if self.func == L.Affine then
    local rel = self.params[1].node_type.value
    assert(rel:isGrid())
    -- TODO: The translated expression for self.params[3] is duplicated.
    assert(self.params[2].m == self.params[2].n + 1)
    assert(self.params[2].n == 3)
    local mat = self.params[2].matvals
    -- Only allowing diagonal translation matrices.
    for i=1,3 do for j=1,3 do
      assert(i == j and mat[i][j] == 1 or
             i ~= j and mat[i][j] == 0)
    end end
    local base = self.params[3]:toRExpr(ctxt)
    local x = mat[1][4]
    local y = mat[2][4]
    local z = mat[3][4]
    if x == 0 and y == 0 and z == 0 then
      return base
    end
    return rexpr (base + {x,y,z}) % [ctxt.relMap[rel]].bounds end
  end
  -- Key unboxing
  -- self.params[1] : AST.Expression
  if self.func == L.id then
    return self.params[1]:toRExpr(ctxt)
  elseif self.func == L.xid then
    return rexpr ( [int3d]([self.params[1]:toRExpr(ctxt)]) ).x end
  elseif self.func == L.yid then
    return rexpr ( [int3d]([self.params[1]:toRExpr(ctxt)]) ).y end
  elseif self.func == L.zid then
    return rexpr ( [int3d]([self.params[1]:toRExpr(ctxt)]) ).z end
  end
  -- Unary arithmetic function
  -- self.params[1] : AST.Expression
  if UNARY_ARITH_FUNS[self.func] then
    local arg = self.params[1]:toRExpr(ctxt)
    return rexpr [UNARY_ARITH_FUNS[self.func]](arg) end
  end
  -- Binary arithmetic function
  -- self.params[1] : AST.Expression
  -- self.params[2] : AST.Expression
  if BINARY_ARITH_FUNS[self.func] then
    local arg1 = self.params[1]:toRExpr(ctxt)
    local arg2 = self.params[2]:toRExpr(ctxt)
    return rexpr [BINARY_ARITH_FUNS[self.func]](arg1, arg2) end
  end
  -- Min/max call
  -- self.params[1] : AST.Expression
  -- self.params[2] : AST.Expression
  if self.func == L.fmax or self.func == L.imax then
    local arg1 = self.params[1]:toRExpr(ctxt)
    local arg2 = self.params[2]:toRExpr(ctxt)
    return rexpr max(arg1, arg2) end
  end
  if self.func == L.fmin or self.func == L.imin then
    local arg1 = self.params[1]:toRExpr(ctxt)
    local arg2 = self.params[2]:toRExpr(ctxt)
    return rexpr min(arg1, arg2) end
  end
  -- Random number generator
  if self.func == L.rand then
    return rexpr [double](C.rand()) / C.RAND_MAX end
  end
  -- Dot product
  -- self.params[1] : AST.Expression
  -- self.params[2] : AST.Expression
  if self.func == L.dot then
    local t1 = self.params[1].node_type
    local t2 = self.params[2].node_type
    assert(t1:isvector() and t2:isvector() and t1.N == t2.N)
    local fun = emitDotProduct(toRType(t1.type), t1.N)
    local arg1 = self.params[1]:toRExpr(ctxt)
    local arg2 = self.params[2]:toRExpr(ctxt)
    return rexpr fun([arg1], [arg2]) end
  end
  -- Element-wise multiplication
  -- self.params[1] : AST.Expression
  -- self.params[2] : AST.Expression
  if self.func == L.times then
    local t1 = self.params[1].node_type
    local t2 = self.params[2].node_type
    assert(t1:isvector() and t2:isvector() and t1.N == t2.N)
    local fun = emitVectorVectorOp('*', toRType(t1.type), t1.N)
    local arg1 = self.params[1]:toRExpr(ctxt)
    local arg2 = self.params[2]:toRExpr(ctxt)
    return rexpr fun([arg1], [arg2]) end
  end
  -- Derive key value from indices
  -- self.params[1] : AST.VectorLiteral
  -- self.params[2] : AST.LuaObject
  --   .node_type.value : R.Relation
  if self.func == L.UNSAFE_ROW then
    local rel = self.params[2].node_type.value
    local vals = newlist(self.params[1].elems)
    vals = vals:map(function(e) return e:toRExpr(ctxt) end)
    if rel:isGrid() then
      assert(#vals == 3)
      return rexpr int3d{ [vals[1]], [vals[2]], [vals[3]] } end
    else
      assert(#vals == 1)
      return rexpr [int1d]([vals[1]]) end
    end
  end
  -- Call to terra function
  -- self.params : AST.Expression*
  if self.func.is_a_terra_func then
    local args = newlist()
    for idx = 1, #self.params do
      args:insert(self.params[idx]:toRExpr(ctxt))
    end
    return rexpr [self.func.terrafn]([args]) end
  end
  -- Print statement
  -- self.params : AST.Expression*
  if self.func == L.print then
    local args = newlist()
    for idx = 1, #self.params do
      local arg = self.params[idx]:toRExpr(ctxt)
      if not self.params[idx].node_type:isvector() then
        args:insert(arg)
      else
        args:insert(rexpr [arg][0] end)
        args:insert(rexpr [arg][1] end)
        args:insert(rexpr [arg][2] end)
      end
    end
    local fmt = ""
    for idx = 1, #self.params do
      local ty = self.params[idx].node_type
      if not ty:isvector() then
        if ty == L.uint64 then
          fmt = fmt .. "%lu "
        elseif ty == L.int then
          fmt = fmt .. "%d "
        else
          fmt = fmt .. "%f "
        end
      else
        local val_fmt
        if ty.type == L.uint64 then
          val_fmt = "%lu"
        elseif ty == L.int then
          fmt = fmt .. "%d"
        else
          val_fmt = "%f"
        end
        fmt = fmt .. "("
        for idx = 1, ty.N do
          fmt = fmt .. val_fmt
          if idx ~= ty.N then
            fmt = fmt .. ", "
          end
        end
        fmt = fmt .. ") "
      end
    end
    fmt = fmt .. "\n"
    return rexpr C.printf([fmt], [args]) end
  end
  -- TODO: Not covered: L.cross, L.length, L.assert
  assert(false)
end
function AST.Cast:toRExpr(ctxt)
  -- self.node_type : T.Type
  -- self.value     : AST.Expression
  return rexpr [toRType(self.node_type)]([self.value:toRExpr(ctxt)]) end
end
function AST.FieldAccess:toRExpr(ctxt)
  -- self.field : R.Field
  -- self.key   : AST.Expression
  local rel = self.key.node_type.relation
  -- TODO: Assuming that off-center accesses are only made using Affine
  -- expressions, therefore we don't need to bounds-check here.
  return rexpr
    [ctxt.relMap[rel]][ [self.key:toRExpr(ctxt)] ].[self.field:Name()]
  end
end
function AST.FieldAccessIndex:toRExpr(ctxt)
  -- self.base  : AST.FieldAccess
  -- self.field : R.Field
  -- self.key   : AST.Expression
  -- self.index : AST.Expression
  return rexpr [self.base:toRExpr(ctxt)][ [self.index:toRExpr(ctxt)] ] end
end
function AST.Global:toRExpr(ctxt)
  -- self.global : PRE.Global
  return rexpr
    [ctxt.globalMap[self.global]]
  end
end
function AST.GlobalIndex:toRExpr(ctxt)
  -- self.index  : AST.Expression
  -- self.global : PRE.Global
  return rexpr [ctxt.globalMap[self.global]][ [self.index:toRExpr(ctxt)] ] end
end
function AST.LetExpr:toRExpr(ctxt)
  -- self.block       : AST.Block
  -- self.exp         : AST.Expression
  -- self.orig_func   : F.Function?
  -- self.orig_params : (AST.Expression*)?
  assert(self.block:is(AST.Block))
  -- Call to user-defined helper function: emit as separate task
  if self.orig_func then
    return recoverHelperCall(self, ctxt)
  end
  -- Call to macro: handle common case of simple alias expression
  assert(#self.block.statements == 1)
  local decl = self.block.statements[1]
  assert(decl:is(AST.DeclStatement))
  ctxt:addAlias(decl.name, decl.initializer:toRExpr(ctxt))
  return self.exp:toRExpr(ctxt)
end
function AST.LuaObject:toRExpr(ctxt)
  quit(self)
end
function AST.MatrixLiteral:toRExpr(ctxt)
  quit(self)
end
function AST.Name:toRExpr(ctxt)
  -- self.name : AST.Symbol
  return assert(ctxt.localMap[self.name])
end
function AST.Number:toRExpr(ctxt)
  -- self.node_type : T.Type
  -- self.value     : number
  return rexpr [toRType(self.node_type)]([self.value]) end
end
function AST.Quote:toRExpr(ctxt)
  -- self.code : AST.Expression
  return self.code:toRExpr(ctxt)
end
function AST.RecordLiteral:toRExpr(ctxt)
  quit(self)
end
function AST.Reduce:toRExpr(ctxt)
  quit(self)
end
function AST.SquareIndex:toRExpr(ctxt)
  -- self.node_type : T.Type
  -- self.base      : AST.Expression
  -- self.index     : AST.Expression
  return rexpr [self.base:toRExpr(ctxt)][ [self.index:toRExpr(ctxt)] ] end
end
function AST.String:toRExpr(ctxt)
  quit(self)
end
function AST.TableLookup:toRExpr(ctxt)
  quit(self)
end
function AST.UnaryOp:toRExpr(ctxt)
  -- self.exp : AST.Expression
  -- self.op  : string
  local t = self.exp.node_type
  local arg = self.exp:toRExpr(ctxt)
  if t:isvector() and self.op == '-' then
    local fun = emitVectorScalarOp('*', toRType(t.type), t.N)
    return rexpr fun(arg, -1) end
  end
  assert(t:isscalar())
  return
    (self.op == '-')   and rexpr -arg    end or
    (self.op == 'not') and rexpr not arg end or
    assert(false)
end
function AST.VectorLiteral:toRExpr(ctxt)
  -- self.node_type : T.Type
  -- self.elems     : AST.Expression*
  return emitVectorCtor(
    toRType(self.node_type.type),
    newlist(self.elems):map(function(e) return e:toRExpr(ctxt) end))
end
function AST.Where:toRExpr(ctxt)
  quit(self)
end

-- FunContext -> RG.rquote
function AST.Statement:toRQuote(ctxt)
  error('Abstract Method')
end
function AST.Assignment:toRQuote(ctxt)
  -- self.lvalue   : AST.Expression
  -- self.exp      : AST.Expression
  -- self.reduceop : string?
  local lval = self.lvalue:toRExpr(ctxt)
  local exp = self.exp:toRExpr(ctxt)
  if self.reduceop then
    return emitReduce(self.reduceop, self.exp.node_type, lval, exp)
  else
    return rquote lval = exp end
  end
end
function AST.Break:toRQuote(ctxt)
  quit(self)
end
function AST.DeclStatement:toRQuote(ctxt)
  -- self.name        : AST.Symbol
  -- self.node_type   : T.Type
  -- self.initializer : AST.Expression
  local rsym = ctxt:addLocal(self.name, self.node_type)
  if self.initializer then
    return rquote var [rsym] = [self.initializer:toRExpr(ctxt)] end
  end
  return rquote var [rsym] end
end
function AST.DeleteStatement:toRQuote(ctxt)
  -- self.key : AST.Expression
  local rel = ctxt.relMap[ctxt.domainRel]
  return rquote
    [rel][ [self.key:toRExpr(ctxt)] ].__valid = false
  end
end
function AST.DoStatement:toRQuote(ctxt)
  -- self.body        : AST.Block
  -- self.orig_func   : F.Function?
  -- self.orig_params : (AST.Expression*)?
  return self.body:toRQuote(ctxt)
end
function AST.ExprStatement:toRQuote(ctxt)
  -- self.exp : AST.Expression | AST.DoStatement
  if self.exp:is(AST.DoStatement) then
    return self.exp:toRQuote(ctxt)
  end
  return rquote [self.exp:toRExpr(ctxt)] end
end
function AST.FieldWrite:toRQuote(ctxt)
  -- self.fieldaccess : AST.FieldAccess
  -- self.exp         : AST.Expression
  -- self.reduceop    : string?
  local lval = self.fieldaccess:toRExpr(ctxt)
  local exp = self.exp:toRExpr(ctxt)
  if self.reduceop then
    return emitReduce(self.reduceop, self.exp.node_type, lval, exp)
  else
    return rquote lval = exp end
  end
end
function AST.GenericFor:toRQuote(ctxt)
  quit(self)
end
function AST.GlobalReduce:toRQuote(ctxt)
  -- self.global   : AST.Global
  -- self.reduceop : string
  -- self.exp      : AST.Expression
  assert(self.global.global == ctxt.reducedGlobal)
  local lval = ctxt.globalReduceAcc
  local exp = self.exp:toRExpr(ctxt)
  return emitReduce(self.reduceop, self.exp.node_type, lval, exp)
end
function AST.IfStatement:toRQuote(ctxt)
  -- self.if_blocks  : AST.CondBlock*
  -- self.else_block : AST.Block
  local quot
  for i=#self.if_blocks,1,-1 do
    local cond = self.if_blocks[i].cond:toRExpr(ctxt)
    local body = self.if_blocks[i].body:toRQuote(ctxt)
    if quot then
      quot = rquote if [cond] then [body] else [quot] end end
    elseif self.else_block then
      local innerElse = self.else_block:toRQuote(ctxt)
      quot = rquote if [cond] then [body] else [innerElse] end end
    else
      quot = rquote if [cond] then [body] end end
    end
  end
  return quot
end
function AST.InsertStatement:toRQuote(ctxt)
  -- self.record   : AST.RecordLiteral
  --   .names      : string*
  --   .exprs      : AST.Expression*
  -- self.relation.node_type.value : R.Relation
  local rg = ctxt.relMap[self.relation.node_type.value]
  local elem = RG.newsymbol(nil, 'elem')
  local fldWriteQuotes = newlist()
  for i,name in ipairs(self.record.names) do
    fldWriteQuotes:insert(rquote
      elem.[name] = [self.record.exprs[i]:toRExpr(ctxt)]
    end)
  end
  return rquote
    var inserted = false
    for [elem] in rg do
      if not elem.__valid then
        elem.__valid = true
        [fldWriteQuotes]
        inserted = true
        break
      end
    end
    RG.assert(inserted, 'Out of space')
  end
end
function AST.NumericFor:toRQuote(ctxt)
  -- self.name  : AST.Symbol
  -- self.lower : AST.Expression
  -- self.upper : AST.Expression
  -- self.body  : AST.Block
  local i = ctxt:addLocal(self.name, L.int)
  return rquote
    for [i] = [self.lower:toRExpr(ctxt)], [self.upper:toRExpr(ctxt)] do
      [self.body:toRQuote(ctxt)]
    end
  end
end
function AST.RepeatStatement:toRQuote(ctxt)
  quit(self)
end
function AST.WhileStatement:toRQuote(ctxt)
  quit(self)
end

-- FunContext -> RG.rquote
function AST.Block:toRQuote(ctxt)
  local stmtQuotes =
    newlist(self.statements):map(function(stmt) return stmt:toRQuote(ctxt) end)
  return rquote
    [stmtQuotes]
  end
end

-------------------------------------------------------------------------------
-- Relation dumping & loading
-------------------------------------------------------------------------------

if USE_HDF then

  local HDF5 = terralib.includec(HDF_HEADER)

  -- HACK: Hardcoding missing #define's
  HDF5.H5F_ACC_TRUNC = 2
  HDF5.H5P_DEFAULT = 0

  -- string*, string, RG.rexpr* -> RG.rquote
  function R.Relation:emitDump(flds, file, vals)
    local flds = flds:copy()
    if self:isCoupled() then
      assert(flds:find(self:CouplingField():Name()))
      flds:insert('__valid')
    end
    local terra create(fname : rawstring, xSize : int, ySize : int, zSize: int)
      var fid = HDF5.H5Fcreate(fname, HDF5.H5F_ACC_TRUNC,
                               HDF5.H5P_DEFAULT, HDF5.H5P_DEFAULT)
      var dataSpace : int32
      escape
        if self:isPlain() or self:isCoupled() then
          emit quote
            var sizes : HDF5.hsize_t[1]
            sizes[0] = xSize
            dataSpace = HDF5.H5Screate_simple(1, sizes, [&uint64](0))
          end
        elseif self:isGrid() then
          emit quote
            -- Legion defaults to column-major layout, so we have to reverse.
            -- This implies that x and z will be flipped in the output file.
            var sizes : HDF5.hsize_t[3]
            sizes[2] = xSize
            sizes[1] = ySize
            sizes[0] = zSize
            dataSpace = HDF5.H5Screate_simple(3, sizes, [&uint64](0))
          end
        else assert(false) end
        local header = newlist() -- terralib.quote*
        local footer = newlist() -- terralib.quote*
        -- terralib.type -> terralib.quote
        local function toHType(T)
          -- TODO: Not supporting: pointers, vectors, non-primitive arrays
          if T:isprimitive() then
            return
              -- HACK: Hardcoding missing #define's
              (T == int)    and HDF5.H5T_STD_I32LE_g  or
              (T == int8)   and HDF5.H5T_STD_I8LE_g   or
              (T == int16)  and HDF5.H5T_STD_I16LE_g  or
              (T == int32)  and HDF5.H5T_STD_I32LE_g  or
              (T == int64)  and HDF5.H5T_STD_I64LE_g  or
              (T == uint)   and HDF5.H5T_STD_U32LE_g  or
              (T == uint8)  and HDF5.H5T_STD_U8LE_g   or
              (T == uint16) and HDF5.H5T_STD_U16LE_g  or
              (T == uint32) and HDF5.H5T_STD_U32LE_g  or
              (T == uint64) and HDF5.H5T_STD_U64LE_g  or
              (T == bool)   and HDF5.H5T_STD_U8LE_g   or
              (T == float)  and HDF5.H5T_IEEE_F32LE_g or
              (T == double) and HDF5.H5T_IEEE_F64LE_g or
              assert(false)
          elseif T:isarray() then
            local elemType = toHType(T.type)
            local arrayType = symbol(HDF5.hid_t, 'arrayType')
            header:insert(quote
              var dims : HDF5.hsize_t[1]
              dims[0] = T.N
              var elemType = [elemType]
              var [arrayType] = HDF5.H5Tarray_create2(elemType, 1, dims)
            end)
            footer:insert(quote
              HDF5.H5Tclose(arrayType)
            end)
            return arrayType
          else assert(false) end
        end
        -- terralib.struct, set(string), string -> ()
        local function emitFieldDecls(fs, whitelist, prefix)
          -- TODO: Only supporting pure structs, not fspaces
          assert(fs:isstruct())
          for _,e in ipairs(fs.entries) do
            local name, type = parseStructEntry(e)
            if whitelist and not whitelist[name] then
              -- do nothing
            elseif type == int2d then
              -- Hardcode special case: int2d structs are stored packed
              local hName = prefix..name
              local int2dType = symbol(HDF5.hid_t, 'int2dType')
              local dataSet = symbol(HDF5.hid_t, 'dataSet')
              header:insert(quote
                var [int2dType] = HDF5.H5Tcreate(HDF5.H5T_COMPOUND, 16)
                HDF5.H5Tinsert(int2dType, "x", 0, HDF5.H5T_STD_I64LE_g)
                HDF5.H5Tinsert(int2dType, "y", 8, HDF5.H5T_STD_I64LE_g)
                var [dataSet] = HDF5.H5Dcreate2(
                  fid, hName, int2dType, dataSpace,
                  HDF5.H5P_DEFAULT, HDF5.H5P_DEFAULT, HDF5.H5P_DEFAULT)
              end)
              footer:insert(quote
                HDF5.H5Dclose(dataSet)
                HDF5.H5Tclose(int2dType)
              end)
            elseif type == int3d then
              -- Hardcode special case: int3d structs are stored packed
              local hName = prefix..name
              local int3dType = symbol(HDF5.hid_t, 'int3dType')
              local dataSet = symbol(HDF5.hid_t, 'dataSet')
              header:insert(quote
                var [int3dType] = HDF5.H5Tcreate(HDF5.H5T_COMPOUND, 24)
                HDF5.H5Tinsert(int3dType, "x", 0, HDF5.H5T_STD_I64LE_g)
                HDF5.H5Tinsert(int3dType, "y", 8, HDF5.H5T_STD_I64LE_g)
                HDF5.H5Tinsert(int3dType, "z", 16, HDF5.H5T_STD_I64LE_g)
                var [dataSet] = HDF5.H5Dcreate2(
                  fid, hName, int3dType, dataSpace,
                  HDF5.H5P_DEFAULT, HDF5.H5P_DEFAULT, HDF5.H5P_DEFAULT)
              end)
              footer:insert(quote
                HDF5.H5Dclose(dataSet)
                HDF5.H5Tclose(int3dType)
              end)
            elseif type:isstruct() then
              emitFieldDecls(type, nil, prefix..name..'.')
            else
              local hName = prefix..name
              local hType = toHType(type)
              local dataSet = symbol(HDF5.hid_t, 'dataSet')
              header:insert(quote
                var hType = [hType]
                var [dataSet] = HDF5.H5Dcreate2(
                  fid, hName, hType, dataSpace,
                  HDF5.H5P_DEFAULT, HDF5.H5P_DEFAULT, HDF5.H5P_DEFAULT)
              end)
              footer:insert(quote
                HDF5.H5Dclose(dataSet)
              end)
            end
          end
        end
        emitFieldDecls(self:fieldSpace(), flds:toSet(), '')
        emit quote [header] end
        emit quote [footer:reverse()] end
      end
      HDF5.H5Sclose(dataSpace)
      HDF5.H5Fclose(fid)
    end
    A.registerFun(create, self:Name()..'_hdf5create_'..flds:join('_'))
    local r = self:regionSymbol()
    local p_r = self:primPartSymbol()
    local s = self:copyRegionSymbol()
    local p_s = self:copyPrimPartSymbol()
    local xSize, ySize, zSize -- RG.rexpr
    if self:isPlain() then
      xSize = rexpr [self:Size():varSymbol()] end
      ySize = rexpr -1 end
      zSize = rexpr -1 end
    elseif self:isGrid() then
      local xNum = self:xNum():varSymbol()
      local yNum = self:yNum():varSymbol()
      local zNum = self:zNum():varSymbol()
      local xBnum = self:xBnum():varSymbol()
      local yBnum = self:yBnum():varSymbol()
      local zBnum = self:zBnum():varSymbol()
      xSize = rexpr xNum + 2*xBnum end
      ySize = rexpr yNum + 2*yBnum end
      zSize = rexpr zNum + 2*zBnum end
    elseif self:isCoupled() then
      xSize = rexpr [self:primPartSize()] * NUM_PRIM_PARTS end
      ySize = rexpr -1 end
      zSize = rexpr -1 end
    else assert(false) end
    return rquote
      var filename = [rawstring](C.malloc(256))
      C.snprintf(filename, 256, file, vals)
      create(filename, xSize, ySize, zSize)
      attach(hdf5, s.[flds], filename, RG.file_read_write)
      for c in [A.primColors()] do
        var p_r_c = p_r[c]
        var p_s_c = p_s[c]
        acquire(p_s_c.[flds])
        copy(p_r_c.[flds], p_s_c.[flds])
        release(p_s_c.[flds])
      end
      detach(hdf5, s.[flds])
      C.free(filename)
    end
  end

  -- string*, string, RG.rexpr* -> RG.rquote
  function R.Relation:emitLoad(flds, file, vals)
    local flds = flds:copy()
    local primPartQuote = rquote end
    if self:isCoupled() then
      assert(flds:find(self:CouplingField():Name()))
      flds:insert('__valid')
      -- XXX: Skipping primary-partition check due to Regent weirdness
      primPartQuote = rquote end
    end
    local r = self:regionSymbol()
    local p_r = self:primPartSymbol()
    local s = self:copyRegionSymbol()
    local p_s = self:copyPrimPartSymbol()
    return rquote
      var filename = [rawstring](C.malloc(256))
      C.snprintf(filename, 256, file, vals)
      attach(hdf5, s.[flds], filename, RG.file_read_only)
      for c in [A.primColors()] do
        var p_r_c = p_r[c]
        var p_s_c = p_s[c]
        acquire(p_s_c.[flds])
        copy(p_s_c.[flds], p_r_c.[flds])
        release(p_s_c.[flds])
      end
      detach(hdf5, s.[flds])
      C.free(filename);
      [primPartQuote]
    end
  end

else -- if not USE_HDF

  -- string*, string, RG.rexpr* -> RG.rquote
  function R.Relation:emitDump(flds, file, vals)
    return rquote
      RG.assert(false, 'HDF I/O not supported; recompile with USE_HDF=1')
    end
  end

  -- string*, string, RG.rexpr* -> RG.rquote
  function R.Relation:emitLoad(flds, file, vals)
    return rquote
      RG.assert(false, 'HDF I/O not supported; recompile with USE_HDF=1')
    end
  end

end -- if USE_HDF

-------------------------------------------------------------------------------
-- Parsing runtime configuration values
-------------------------------------------------------------------------------

local JSON = terralib.includec('json.h')

local Enum = {}

-- string* -> Enum
function A.Enum(...)
  local enum = setmetatable({}, Enum)
  enum.__choices = {...}
  for i,val in ipairs({...}) do
    assert(type(val) == 'string')
    enum[val] = i-1
  end
  return enum
end

-- A -> bool
local function equalsDoubleVec(x)
  return terralib.types.istype(x) and x:isarray() and x.type == double
end

-- ConfigMap = map(string,ConfigKind)
-- ConfigKind = Enum | int | double | double[N] | ConfigMap

-- () -> RG.symbol
A.configSymbol = terralib.memoize(function()
  return RG.newsymbol(nil, 'config')
end)

-- ConfigMap
local CONFIG_MAP = {}

-- () -> terralib.struct
A.configStruct = terralib.memoize(function()
  local function fillStruct(struct_, map)
    for fld,kind in pairs(map) do
      if getmetatable(kind) == Enum then
        struct_.entries:insert({field=fld, type=int})
      elseif kind == int or kind == double or equalsDoubleVec(kind) then
        struct_.entries:insert({field=fld, type=kind})
      else
        local nested = terralib.types.newstruct(fld)
        fillStruct(nested, kind)
        struct_.entries:insert({field=fld, type=nested})
      end
    end
    A.registerStruct(struct_)
  end
  local s = terralib.types.newstruct('Config')
  fillStruct(s, CONFIG_MAP)
  return s
end)

-- string, terralib.expr? -> terralib.quote
local function errorOut(msg, name)
  if name then
    return quote
      C.printf('%s for option %s\n', msg, name)
      C.exit(1)
    end
  else
    return quote
      C.printf('%s\n', msg)
      C.exit(1)
    end
  end
end

-- terralib.symbol, terralib.expr, terralib.expr, ConfigKind -> terralib.quote
local function emitValueParser(name, lval, rval, kind)
  if getmetatable(kind) == Enum then
    return quote
      if [rval].type ~= JSON.json_string then
        [errorOut('Wrong type', name)]
      end
      var found = false
      escape for i,choice in ipairs(kind.__choices) do emit quote
        if C.strcmp([rval].u.string.ptr, choice) == 0 then
          [lval] = i-1
          found = true
        end
      end end end
      if not found then
        [errorOut('Unexpected value', name)]
      end
    end
  elseif kind == int then
    return quote
      if [rval].type ~= JSON.json_integer then
        [errorOut('Wrong type', name)]
      end
      [lval] = [rval].u.integer
    end
  elseif kind == double then
    return quote
      if [rval].type ~= JSON.json_double then
        [errorOut('Wrong type', name)]
      end
      [lval] = [rval].u.dbl
    end
  elseif equalsDoubleVec(kind) then
    return quote
      if [rval].type ~= JSON.json_array then
        [errorOut('Wrong type', name)]
      end
      if [rval].u.array.length ~= [kind.N] then
        [errorOut('Wrong length', name)]
      end
      for i = 0,[kind.N] do
        var rval_i = [rval].u.array.values[i]
        if rval_i.type ~= JSON.json_double then
          [errorOut('Wrong element type', name)]
        end
        [lval][i] = rval_i.u.dbl
      end
    end
  else
    return quote
      var totalParsed = 0
      if [rval].type ~= JSON.json_object then
        [errorOut('Wrong type', name)]
      end
      for i = 0,[rval].u.object.length do
        var nodeName = [rval].u.object.values[i].name
        var nodeValue = [rval].u.object.values[i].value
        var parsed = false
        escape for fld,subKind in pairs(kind) do emit quote
          if C.strcmp(nodeName, fld) == 0 then
            [emitValueParser(nodeName, `[lval].[fld], nodeValue, subKind)]
            parsed = true
          end
        end end end
        if parsed then
          totalParsed = totalParsed + 1
        else
          C.printf('Ignoring option %s\n', nodeName)
        end
      end
      -- TODO: Assuming the json file contains no duplicate values
      if totalParsed < [UTIL.tableSize(kind)] then
        [errorOut('Missing options from config file')]
      end
    end
  end
end

-- () -> (&int8 -> ...)
local emitConfigParser = terralib.memoize(function()
  local configStruct = A.configStruct()
  local terra parseConfig(filename : &int8) : configStruct
    var config : configStruct
    var f = C.fopen(filename, 'r');
    if f == nil then [errorOut('Cannot open config file')] end
    var res1 = C.fseek(f, 0, C.SEEK_END);
    if res1 ~= 0 then [errorOut('Cannot seek to end of config file')] end
    var len = C.ftell(f);
    if len < 0 then [errorOut('Cannot ftell config file')] end
    var res2 = C.fseek(f, 0, C.SEEK_SET);
    if res2 ~= 0 then [errorOut('Cannot seek to start of config file')] end
    var buf = [&int8](C.malloc(len));
    if buf == nil then [errorOut('Malloc error while parsing config')] end
    var res3 = C.fread(buf, 1, len, f);
    if res3 < len then [errorOut('Cannot read from config file')] end
    C.fclose(f);
    var errMsg : int8[256]
    var settings = JSON.json_settings{ 0, 0, nil, nil, nil, 0 }
    settings.settings = JSON.json_enable_comments
    var root = JSON.json_parse_ex(&settings, buf, len, errMsg)
    if root == nil then
      C.printf('JSON parsing error: %s\n', errMsg)
      C.exit(1)
    end
    [emitValueParser('<root>', config, root, CONFIG_MAP)]
    JSON.json_value_free(root)
    C.free(buf)
    return config
  end
  A.registerFun(parseConfig, 'parseConfig')
  return parseConfig
end)

-- string, Enum | int | double | double[N] -> M.AST.Expr
function A.readConfig(name, type)
  assert(getmetatable(type) == Enum or type == int or type == double or
         equalsDoubleVec(type))
  local map = CONFIG_MAP
  local fields = name:split('.')
  for i,fld in ipairs(fields) do
    if i == #fields then
      if map[fld] then assert(map[fld] == type) else map[fld] = type end
    else
      if not map[fld] then map[fld] = {} end
      map = map[fld]
    end
  end
  return M.AST.ReadConfig(name)
end

-- string, Enum | int | double | double[N] -> PRE.Global
function A.globalFromConfig(name, type)
  local lType =
    getmetatable(type) == Enum and L.int                      or
    type == int                and L.int                      or
    type == double             and L.double                   or
    equalsDoubleVec(type)      and L.vector(L.double, type.N) or
    assert(false)
  return L.Global(name, lType, A.readConfig(name, type))
end

-------------------------------------------------------------------------------
-- Control program translation
-------------------------------------------------------------------------------

-- () -> RG.rquote
function M.AST.Stmt:toRQuote()
  error('Abstract method')
end
function M.AST.Block:toRQuote()
  return rquote
    [self.stmts:map(function(s) return s:toRQuote() end)]
  end
end
function M.AST.ForEach:toRQuote()
  assert(not self.subset)
  -- Translate kernel to task
  local info = self.fun:_get_typechecked(42, self.rel, {})
  local tsk, fCtxt = self.fun:toKernelTask()
  -- Collect arguments for call
  local actualArgs = newlist()
  local c = RG.newsymbol(nil, 'c')
  if mustParallelizeManually(info) then
    for _,rel in ipairs(fCtxt.accessedRels) do
      actualArgs:insert(rexpr [rel:primPartSymbol()][c] end)
    end
  else
    for _,rel in ipairs(fCtxt.accessedRels) do
      actualArgs:insert(rel:regionSymbol())
    end
  end
  for _,g in ipairs(fCtxt.readGlobals) do
    actualArgs:insert(g:varSymbol())
  end
  -- Synthesize call expression
  local callExpr = rexpr [tsk]([actualArgs]) end
  local callQuote = rquote [callExpr] end
  if fCtxt.reducedGlobal then
    local retSym = fCtxt.reducedGlobal:varSymbol()
    local op = fCtxt.globalReduceOp
    callQuote =
      (op == '+')   and rquote [retSym] +=   [callExpr]     end or
      (op == '-')   and rquote [retSym] +=   [callExpr]     end or
      (op == '*')   and rquote [retSym] *=   [callExpr]     end or
      (op == '/')   and rquote [retSym] *=   [callExpr]     end or
      (op == 'max') and rquote [retSym] max= [callExpr]     end or
      (op == 'min') and rquote [retSym] min= [callExpr]     end or
      assert(false)
  end
  if mustParallelizeManually(info) then
    callQuote = rquote
      for [c] in [A.primColors()] do
        [callQuote]
      end
    end
  end
  -- Update primary partitioning (if appropriate)
  local primPartQuote = rquote end
  if self.rel:isCoupled() then
    local pt = info.field_use[self.rel:CouplingField()]
    if pt and pt.write then
      primPartQuote = self.rel:emitPrimPartUpdate()
    end
  end
  return rquote
    [callQuote];
    [primPartQuote]
  end
end
function M.AST.If:toRQuote()
  if self.elseBlock then
    return rquote
      if [self.cond:toRExpr()] then
        [self.thenBlock:toRQuote()]
      else
        [self.elseBlock:toRQuote()]
      end
    end
  else
    return rquote
      if [self.cond:toRExpr()] then
        [self.thenBlock:toRQuote()]
      end
    end
  end
end
function M.AST.FillField:toRQuote()
  error('Fill operations are handled separately')
end
function M.AST.SetGlobal:toRQuote()
  return rquote
    [self.global:varSymbol()] = [self.expr:toRExpr()]
  end
end
function M.AST.While:toRQuote()
  if self.spmd then
    return rquote
      __demand(__spmd)
      while [self.cond:toRExpr()] do
        [self.body:toRQuote()]
      end
    end
  else
    return rquote
      while [self.cond:toRExpr()] do
        [self.body:toRQuote()]
      end
    end
  end
end
function M.AST.Do:toRQuote()
  if self.spmd then
    return rquote
      __demand(__spmd)
      do
        [self.body:toRQuote()]
      end
    end
  else
    return rquote
      do
        [self.body:toRQuote()]
      end
    end
  end
end
function M.AST.Print:toRQuote()
  local formals =
    self.globals:map(function(g) return RG.newsymbol(toRType(g:Type())) end)
  local task output([formals])
    C.printf([self.fmt], [formals])
  end
  A.registerTask(output, 'output')
  local actuals = self.globals:map(function(g) return g:varSymbol() end)
  return rquote
    output([actuals])
  end
end
function M.AST.Dump:toRQuote()
  local valRExprs = self.vals:map(function(v) return v:toRExpr() end)
  return self.rel:emitDump(self.flds, self.file, valRExprs)
end
function M.AST.Load:toRQuote()
  local valRExprs = self.vals:map(function(v) return v:toRExpr() end)
  return self.rel:emitLoad(self.flds, self.file, valRExprs)
end
function M.AST.Inline:toRQuote()
  return self.quot
end
function M.AST.Error:toRQuote()
  return rquote RG.assert(false, [self.msg]) end
end

-- () -> RG.rexpr
function M.AST.Cond:toRExpr()
  error('Abstract method')
end
function M.AST.Literal:toRExpr()
  return rexpr [self.val] end
end
function M.AST.And:toRExpr()
  return rexpr [self.lhs:toRExpr()] and [self.rhs:toRExpr()] end
end
function M.AST.Or:toRExpr()
  return rexpr [self.lhs:toRExpr()] or [self.rhs:toRExpr()] end
end
function M.AST.Not:toRExpr()
  return rexpr not [self.cond:toRExpr()] end
end
function M.AST.Compare:toRExpr()
  local a = self.lhs:toRExpr()
  local b = self.rhs:toRExpr()
  return
    (self.op == '==') and rexpr a == b end or
    (self.op == '~=') and rexpr a ~= b end or
    (self.op == '>')  and rexpr a > b  end or
    (self.op == '<')  and rexpr a < b  end or
    (self.op == '>=') and rexpr a >= b end or
    (self.op == '<=') and rexpr a <= b end or
    assert(false)
end

-- () -> RG.rexpr
function M.AST.Expr:toRExpr()
  error('Abstract method')
end
function M.AST.Const:toRExpr()
  return rexpr [toRConst(self.val)] end
end
function M.AST.GetGlobal:toRExpr()
  return rexpr [self.global:varSymbol()] end
end
function M.AST.BinaryOp:toRExpr()
  local a = self.lhs:toRExpr()
  local b = self.rhs:toRExpr()
  return
    (self.op == '+')   and rexpr    a + b end or
    (self.op == '-')   and rexpr    a - b end or
    (self.op == '*')   and rexpr    a * b end or
    (self.op == '/')   and rexpr    a / b end or
    (self.op == '%')   and rexpr    a % b end or
    (self.op == 'max') and rexpr max(a,b) end or
    (self.op == 'min') and rexpr min(a,b) end or
    assert(false)
end
function M.AST.UnaryOp:toRExpr()
  local a = self.arg:toRExpr()
  return
    (self.op == '-') and rexpr -a end or
    assert(false)
end
function M.AST.Array:toRExpr()
  return rexpr array([self.elems:map(function(e) return e:toRExpr() end)]) end
end
function M.AST.Index:toRExpr()
  return rexpr [self.base:toRExpr()][ [self.index] ] end
end
function M.AST.ReadConfig:toRExpr()
  local expr = rexpr [A.configSymbol()] end
  local fields = self.name:split('.')
  for _,fld in ipairs(fields) do
    expr = rexpr expr.[fld] end
  end
  return expr
end
function M.AST.Cond2Expr:toRExpr()
  return self.cond:toRExpr()
end

-- M.AST.FillField -> RG.rquote*
local function emitFillTaskCalls(fillStmts)
  local fillsPerRelation = {} -- map(R.Relation,M.AST.FillField)
  for _,fill in ipairs(fillStmts) do
    local rel = fill.fld:Relation()
    local stmts = fillsPerRelation[rel] or newlist()
    stmts:insert(fill)
    fillsPerRelation[rel] = stmts
    assert(not (rel:isCoupled() and fill.fld == rel:CouplingField()))
  end
  local fillTasks = {} -- map(R.Relation,RG.task)
  for rel,fills in pairs(fillsPerRelation) do
    local arg = RG.newsymbol(rel:regionType(), rel:Name())
    local privileges = newlist()
    privileges:insert(RG.privilege(RG.reads, arg))
    privileges:insert(RG.privilege(RG.writes, arg))
    local body
    if not rel:isCoupled() and RG.config['openmp'] then
      body = fills:map(function(fill) return rquote
        __demand(__openmp)
        for e in arg do
          e.[fill.fld:Name()] = [toRConst(fill.val, fill.fld:Type())]
        end
      end end)
    else
      body = fills:map(function(fill) return rquote
        for e in arg do
          e.[fill.fld:Name()] = [toRConst(fill.val, fill.fld:Type())]
        end
      end end)
    end
    local tsk
    if RG.check_cuda_available() then
      __demand(__parallel, __cuda)
      task tsk([arg]) where [privileges] do
        [body]
      end
    else
      __demand(__parallel)
      task tsk([arg]) where [privileges] do
        [body]
      end
    end
    fillTasks[rel] = tsk
    A.registerTask(tsk, rel:Name()..'_fillTask')
  end
  local calls = newlist() -- RG.rquote*
  for rel,tsk in pairs(fillTasks) do
    calls:insert(rquote [tsk]([rel:regionSymbol()]) end)
  end
  return calls
end

-- () -> RG.symbol()
PRE.Global.varSymbol = terralib.memoize(function(self)
  return RG.newsymbol(toRType(self:Type()), idSanitize(self:Name()))
end)

-- M.AST.Expr, M.AST.Expr, M.AST.Expr -> ()
function A.translate(xTiles, yTiles, zTiles)
  local header = newlist() -- RG.rquote*
  local body = newlist() -- RG.rquote*
  -- Collect declarations
  local globals = newlist() -- PRE.Global*
  local globalInits = {} -- map(PRE.Global, M.AST.Expr)
  local rels = newlist() -- R.Relation*
  for _,decl in ipairs(M.decls()) do
    if M.AST.NewField.check(decl) then
      -- Do nothing
    elseif M.AST.NewFunction.check(decl) then
      -- Do nothing
    elseif M.AST.NewGlobal.check(decl) then
      globals:insert(decl.global)
      globalInits[decl.global] = decl.init
    elseif M.AST.NewRelation.check(decl) then
      rels:insert(decl.rel)
    else assert(false) end
  end
  -- Emit global declarations
  header:insert(rquote
    var [NX] = [xTiles:toRExpr()]
    var [NY] = [yTiles:toRExpr()]
    var [NZ] = [zTiles:toRExpr()]
  end)
  for _,g in ipairs(globals) do
    header:insert(rquote
      var [g:varSymbol()] = [globalInits[g]:toRExpr()]
    end)
  end
  -- Emit region declarations
  for _,rel in ipairs(rels) do
    header:insert(rel:emitRegionInit())
  end
  -- Emit primary partitioning scheme
  header:insert(rquote
    var [A.primColors()] = ispace(int3d, {NX,NY,NZ})
  end)
  for _,rel in ipairs(rels) do
    header:insert(rel:emitPrimPartInit())
    -- Emit transfer queues for coupled relations
    if rel:isCoupled() then
      for i = 1,#rel:XferStencil() do
        header:insert(rel:emitQueueInit(i))
        header:insert(rel:emitQueuePartInit(i))
      end
    end
  end
  -- Process fill statements
  local fillStmts = newlist()
  for _,s in ipairs(M.stmts()) do
    if M.AST.FillField.check(s) then
      fillStmts:insert(s)
    end
  end
  body:insertall(emitFillTaskCalls(fillStmts))
  for _,rel in ipairs(rels) do
    if rel:isCoupled() then
      body:insert(rel:emitValidInit())
    end
  end
  -- Process other statements
  for _,s in ipairs(M.stmts()) do
    if not M.AST.FillField.check(s) then
      body:insert(s:toRQuote())
    end
  end
  -- Synthesize work task
  local opts = newlist() -- RG.rexpr*
  opts:insertall(rels:map(function(rel) return rel:primPartSymbol() end))
  opts:insert(A.primColors())
  for _,domRel in ipairs(rels) do
    local function recordSubsumption(connFld)
      local rngRel = connFld:Type().relation
      local domRg = domRel:regionSymbol()
      local domPart = domRel:primPartSymbol()
      local rngRg = rngRel:regionSymbol()
      local rngPart = rngRel:primPartSymbol()
      opts:insert(rexpr
        image(rngRg, domPart, domRg.[connFld:Name()]) <= rngPart
      end)
    end
    if domRel:isCoupled() then
      recordSubsumption(domRel:CouplingField())
    end
    if domRel:isGrid() then
      for _,fld in ipairs(domRel:CoarseningFields()) do
        recordSubsumption(fld)
      end
    end
  end
  A.configSymbol():settype(A.configStruct())
  local task work([A.configSymbol()])
    [header]
    __parallelize_with [opts] do [body] end
  end
  A.registerTask(work, 'work')
  -- Synthesize main task
  local parseConfig = emitConfigParser()
  local task main()
    var args = RG.c.legion_runtime_get_input_args()
    for i = 1,args.argc do
      if C.strcmp(args.argv[i], '-i') == 0 and i < args.argc - 1 then
        work(parseConfig(args.argv[i+1]))
      end
    end
  end
  A.registerTask(main, 'main')
  -- Emit to executable
  RG.saveobj(main, OBJNAME, 'executable', nil, LIBS)
end
