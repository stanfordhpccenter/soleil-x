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

local R = {}
package.loaded["ebb.src.relations"] = R

local Pre   = require "ebb.src.prelude"
local T     = require "ebb.src.types"
local F     = require "ebb.src.functions"
local M     = require "ebb.src.main"
local UTIL  = require "ebb.src.util"

local all        = UTIL.all

-------------------------------------------------------------------------------

local valid_name_err_msg_base =
  "must be valid Lua Identifiers: a letter or underscore,"..
  " followed by zero or more underscores, letters, and/or numbers"
local valid_name_err_msg = {
  relation = "Relation names "  .. valid_name_err_msg_base,
  field    = "Field names "     .. valid_name_err_msg_base,
  subset   = "Subset names "    .. valid_name_err_msg_base
}
local function is_valid_lua_identifier(name)
  if type(name) ~= 'string' then return false end
  -- regex for valid LUA identifiers
  if not name:match('^[_%a][_%w]*$') then return false end
  return true
end

local function is_num(x) return type(x) == 'number' end

local function is_bool(x) return type(x) == 'boolean' end

local function is_int(x) return type(x) == 'number' and x == math.floor(x) end

local function is_int_global(x)
  return Pre.is_global(x) and x:Type() == T.int
end

local function is_double_global(x)
  return Pre.is_global(x) and x:Type() == T.double
end

local function is_bool_global(x)
  return Pre.is_global(x) and x:Type() == T.bool
end

-------------------------------------------------------------------------------

local Relation    = {}
Relation.__index  = Relation
R.Relation        = Relation
local function is_relation(obj) return getmetatable(obj) == Relation end
R.is_relation     = is_relation


local Field       = {}
Field.__index     = Field
R.Field           = Field
local function is_field(obj) return getmetatable(obj) == Field end
R.is_field        = is_field

local Subset      = {}
Subset.__index    = Subset
R.Subset          = Subset
local function is_subset(obj) return getmetatable(obj) == Subset end
R.is_subset       = is_subset

-------------------------------------------------------------------------------
--[[  Relation methods                                                     ]]--
-------------------------------------------------------------------------------

-- A Relation must be in one of the following MODES:
--    PLAIN
--    GRID
--    COUPLED
function Relation:isPlain()       return self._mode == 'PLAIN'      end
function Relation:isGrid()        return self._mode == 'GRID'       end
function Relation:isCoupled()     return self._mode == 'COUPLED'    end

local errorMsg = [[NewRelation must be called with the following parameters:
  name : string
  mode : 'PLAIN' | 'GRID' | 'COUPLED'
  -- if mode == 'PLAIN':
  size : L.Global(L.int)
  -- if mode == 'GRID':
  xNum, yNum, zNum : L.Global(L.int)
  xOrigin, yOrigin, zOrigin : L.Global(L.double)
  xWidth, yWidth, zWidth : L.Global(L.double)
  xBnum, yBnum, zBnum : L.Global(L.int)
  xPeriodic, yPeriodic, zPeriodic : L.Global(L.bool)
  -- if mode == 'COUPLED':
  size : L.Global(L.int)
  coupled_with : Relation
  coupling_field : string
  max_skew : L.Global(L.double)
  max_xfer_num : L.Global(L.int)
  xfer_stencil : int[3]*]]

local relation_uid = 0
function R.NewRelation(params)
  -- Check the parameters
  local function checkParams(params)
    local check = true
      and type(params) == 'table'
      and is_valid_lua_identifier(params.name)
      and (params.mode == 'PLAIN' or
           params.mode == 'GRID' or
           params.mode == 'COUPLED')
    if params.mode == 'PLAIN' then
      check = check
        and params.size and is_int_global(params.size)
    elseif params.mode == 'GRID' then
      check = check
        and params.xNum and is_int_global(params.xNum)
        and params.yNum and is_int_global(params.yNum)
        and params.zNum and is_int_global(params.zNum)
        and params.xOrigin and is_double_global(params.xOrigin)
        and params.yOrigin and is_double_global(params.yOrigin)
        and params.zOrigin and is_double_global(params.zOrigin)
        and params.xWidth and is_double_global(params.xWidth)
        and params.yWidth and is_double_global(params.yWidth)
        and params.zWidth and is_double_global(params.zWidth)
        and params.xBnum and is_int_global(params.xBnum)
        and params.yBnum and is_int_global(params.yBnum)
        and params.zBnum and is_int_global(params.zBnum)
        and params.xPeriodic and is_bool_global(params.xPeriodic)
        and params.yPeriodic and is_bool_global(params.yPeriodic)
        and params.zPeriodic and is_bool_global(params.zPeriodic)
    elseif params.mode == 'COUPLED' then
      check = check
        and params.size and is_int_global(params.size)
        and is_relation(params.coupled_with)
        and is_valid_lua_identifier(params.coupling_field)
        and params.max_skew and is_double_global(params.max_skew)
        and params.max_xfer_num and is_int_global(params.max_xfer_num)
        and terralib.israwlist(params.xfer_stencil)
        and all(params.xfer_stencil, function(dir) return
                  terralib.israwlist(dir) and #dir == 3 and all(dir, is_int)
                end)
    else assert(false) end
    return check
  end
  if not checkParams(params) then error(errorMsg, 2) end

  -- Construct the relation
  local rel = setmetatable( {
    _name       = params.name,
    _mode       = params.mode,
    _uid        = relation_uid,
    _fields     = terralib.newlist(),
    _macros     = terralib.newlist(),
    _functions  = terralib.newlist(),
  }, Relation)
  relation_uid = relation_uid + 1 -- increment unique id counter

  -- Perform mode-dependent setup
  if params.mode == 'PLAIN' then
    rawset(rel, '_size',  params.size)
  elseif params.mode == 'GRID' then
    rawset(rel, '_xNum', params.xNum)
    rawset(rel, '_yNum', params.yNum)
    rawset(rel, '_zNum', params.zNum)
    rawset(rel, '_xOrigin', params.xOrigin)
    rawset(rel, '_yOrigin', params.yOrigin)
    rawset(rel, '_zOrigin', params.zOrigin)
    rawset(rel, '_xWidth', params.xWidth)
    rawset(rel, '_yWidth', params.yWidth)
    rawset(rel, '_zWidth', params.zWidth)
    rawset(rel, '_xBnum', params.xBnum)
    rawset(rel, '_yBnum', params.yBnum)
    rawset(rel, '_zBnum', params.zBnum)
    rawset(rel, '_xPeriodic', params.xPeriodic)
    rawset(rel, '_yPeriodic', params.yPeriodic)
    rawset(rel, '_zPeriodic', params.zPeriodic)
    rawset(rel, '_coarsening_fields', terralib.newlist())
  elseif params.mode == 'COUPLED' then
    rawset(rel, '_size',  params.size)
    rawset(rel, '_coupling_field',
           rel:NewField(params.coupling_field, params.coupled_with))
    rawset(rel, '_max_skew', params.max_skew)
    rawset(rel, '_max_xfer_num', params.max_xfer_num)
    rawset(rel, '_xfer_stencil', terralib.newlist(params.xfer_stencil))
  else assert(false) end

  -- Register & return the relation
  M.decls():insert(M.AST.NewRelation(rel))
  return rel
end

function Relation:_INTERNAL_UID()
  return self._uid
end

function Relation:Name()
  return self._name
end

function Relation:Mode()
  return self._mode
end

function Relation:Size()
  assert(self:isPlain() or self:isCoupled())
  return self._size
end

function Relation:NumDims()
  return
    self:isPlain()   and 1 or
    self:isGrid()    and 3 or
    self:isCoupled() and 1 or
    assert(false)
end

function Relation:Fields()
  return self._fields
end

function Relation:foreach(user_func, ...)
  if not F.is_function(user_func) then
    error('foreach(): expects an ebb function as the first argument', 2)
  end
  user_func:_doForEach(self, ...)
end

-- prevent user from modifying the lua table
function Relation:__newindex(fieldname,value)
  error("Cannot assign members to Relation object "..
      "(did you mean to call relation:New...?)", 2)
end

-------------------------------------------------------------------------------
--[[  Grids                                                                ]]--
-------------------------------------------------------------------------------

function R.Relation:xNum()
  assert(self:isGrid())
  return self._xNum
end
function R.Relation:yNum()
  assert(self:isGrid())
  return self._yNum
end
function R.Relation:zNum()
  assert(self:isGrid())
  return self._zNum
end

function R.Relation:xOrigin()
  assert(self:isGrid())
  return self._xOrigin
end
function R.Relation:yOrigin()
  assert(self:isGrid())
  return self._yOrigin
end
function R.Relation:zOrigin()
  assert(self:isGrid())
  return self._zOrigin
end

function R.Relation:xWidth()
  assert(self:isGrid())
  return self._xWidth
end
function R.Relation:yWidth()
  assert(self:isGrid())
  return self._yWidth
end
function R.Relation:zWidth()
  assert(self:isGrid())
  return self._zWidth
end

function R.Relation:xBnum()
  assert(self:isGrid())
  return self._xBnum
end
function R.Relation:yBnum()
  assert(self:isGrid())
  return self._yBnum
end
function R.Relation:zBnum()
  assert(self:isGrid())
  return self._zBnum
end

function R.Relation:xPeriodic()
  assert(self:isGrid())
  return self._xPeriodic
end
function R.Relation:yPeriodic()
  assert(self:isGrid())
  return self._yPeriodic
end
function R.Relation:zPeriodic()
  assert(self:isGrid())
  return self._zPeriodic
end

function R.Relation:CoarseningFields()
  assert(self:isGrid())
  return self._coarsening_fields:copy()
end

-------------------------------------------------------------------------------
--[[  Coupled relations                                                    ]]--
-------------------------------------------------------------------------------

function Relation:CouplingField()
  assert(self:isCoupled())
  return self._coupling_field
end

function Relation:MaxSkew()
  assert(self:isCoupled())
  return self._max_skew
end

function Relation:MaxXferNum()
  assert(self:isCoupled())
  return self._max_xfer_num
end

function Relation:XferStencil()
  assert(self:isCoupled())
  return self._xfer_stencil
end

-------------------------------------------------------------------------------
--[[  Virtual fields                                                       ]]--
-------------------------------------------------------------------------------

function Relation:NewFieldMacro(name, macro)
  if not name or type(name) ~= "string" then
    error("NewFieldMacro() expects a string as the first argument", 2)
  end
  if not is_valid_lua_identifier(name) then
    error(valid_name_err_msg.field, 2)
  end
  if self[name] then
    error("Cannot create a new field-macro with name '"..name.."'  "..
          "That name is already being used.", 2)
  end

  if not Pre.is_macro(macro) then
    error("NewFieldMacro() expects a Macro as the 2nd argument", 2)
  end

  rawset(self, name, macro)
  self._macros:insert(macro)
  return macro
end

local FieldDispatcher     = {}
FieldDispatcher.__index   = FieldDispatcher
R.FieldDispatcher         = FieldDispatcher
local function NewFieldDispatcher()
  return setmetatable({
    _reader   = nil,
    _writer   = nil,
    _reducers = {},
  }, FieldDispatcher)
end
local function isfielddispatcher(obj)
  return getmetatable(obj) == FieldDispatcher
end
R.isfielddispatcher = isfielddispatcher

local function getFieldDispatcher(rel, fname, ufunc)
  if not fname or type(fname) ~= "string" then
    error("NewField*Function() expects a string as the first argument", 3)
  end
  if not is_valid_lua_identifier(fname) then
    error(valid_name_err_msg.field, 3)
  end
  if not F.is_function(ufunc) then
    error("NewField*Function() expects an Ebb Function "..
          "as the last argument", 3)
  end

  local lookup = rel[fname]
  if lookup and isfielddispatcher(lookup) then return lookup
  elseif lookup then
    error("Cannot create a new field-function with name '"..fname.."'  "..
          "That name is already being used.", 3)
  end

  rawset(rel, fname, NewFieldDispatcher())
  return rel[fname]
end

function Relation:NewFieldReadFunction(fname, userfunc)
  local dispatch = getFieldDispatcher(self, fname, userfunc)
  if dispatch._reader then
    error("NewFieldReadFunction() error: function already assigned.", 2)
  end
  dispatch._reader = userfunc
  self._functions:insert(userfunc)
  return userfunc
end

function Relation:NewFieldWriteFunction(fname, userfunc)
  local dispatch = getFieldDispatcher(self, fname, userfunc)
  if dispatch._writer then
    error("NewFieldWriteFunction() error: function already assigned.", 2)
  end
  dispatch._writer = userfunc
  self._functions:insert(userfunc)
  return userfunc
end

local redops = {
  ['+'] = true,
  ['-'] = true,
  ['*'] = true,
  ['max'] = true,
  ['min'] = true,
}
function Relation:NewFieldReduceFunction(fname, op, userfunc)
  local dispatch = getFieldDispatcher(self, fname, userfunc)
  if not redops[op] then
    error("NewFieldReduceFunction() expects a reduction operator as the "..
          "second argument.", 2)
  end
  if dispatch._reducers[op] then
    error("NewFieldReduceFunction() error: '"..op.."' "..
          "function already assigned.", 2)
  end
  dispatch._reducers[op] = userfunc
  self._functions:insert(userfunc)
  return userfunc
end

-------------------------------------------------------------------------------
--[[  Subsets                                                              ]]--
-------------------------------------------------------------------------------

function Subset:foreach(user_func, ...)
  if not F.is_function(user_func) then
    error('map(): expects an Ebb function as the argument', 2)
  end
  user_func:_doForEach(self, ...)
end

function Subset:Relation()
  return self._owner
end

function Subset:Name()
  return self._name
end

function Subset:FullName()
  return self._owner._name .. '.' .. self._name
end

function Subset:Rectangle()
  return self._rectangle
end

-- prevent user from modifying the lua table
function Subset:__newindex(name,value)
  error("Cannot assign members to Subset object", 2)
end

local function is_int(obj)
  return type(obj) == 'number' and obj % 1 == 0
end

-------------------------------------------------------------------------------
--[[  Fields                                                               ]]--
-------------------------------------------------------------------------------

-- prevent user from modifying the lua table
function Field:__newindex(name,value)
  error("Cannot assign members to Field object", 2)
end

function Field:Name()
  return self._name
end
function Field:FullName()
  return self._owner._name .. '.' .. self._name
end
function Field:Size()
  return self._owner:Size()
end
function Field:Type()
  return self._type
end
function Field:Relation()
  return self._owner
end

function Relation:NewField(name, typ)
  if not name or type(name) ~= "string" then
    error("NewField() expects a string as the first argument", 2)
  end
  if not is_valid_lua_identifier(name) then
    error(valid_name_err_msg.field, 2)
  end
  if self[name] then
    error("Cannot create a new field with name '"..name.."'  "..
          "That name is already being used.", 2)
  end
  if string.sub(name, 1, 2) == '__' then
    error("Field names starting with '__' are reserved by the system.", 2)
  end

  if is_relation(typ) then
    typ = T.key(typ)
  end
  if not T.istype(typ) or not typ:isfieldvalue() then
    error("NewField() expects an Ebb type or "..
          "relation as the 2nd argument", 2)
  end
  if typ:iskey() and typ.relation:isCoupled() then
    error("Foreign keys to coupled relations are not allowed", 2)
  end

  -- create the field
  local field   = setmetatable({
    _type   = typ,
    _name   = name,
    _owner  = self,
  }, Field)
  rawset(self, name, field)

  self._fields:insert(field)

  M.decls():insert(M.AST.NewField(field))

  return field
end

function Field:Fill(val)
  if M.isExprConst(val) then val = M.AST.Const(val) end
  M.stmts():insert(M.AST.FillField(self, val))
end

-------------------------------------------------------------------------------
--[[  Input/output                                                         ]]--
-------------------------------------------------------------------------------

function Relation:Dump(flds, file, ...)
  local args = terralib.newlist({...}):map(function(x)
    return M.isExprConst(x) and M.AST.Const(x) or x
  end)
  M.stmts():insert(M.AST.Dump(self, terralib.newlist(flds), file, args))
end

function Relation:Load(flds, file, ...)
  local args = terralib.newlist({...}):map(function(x)
    return M.isExprConst(x) and M.AST.Const(x) or x
  end)
  M.stmts():insert(M.AST.Load(self, terralib.newlist(flds), file, args))
end
