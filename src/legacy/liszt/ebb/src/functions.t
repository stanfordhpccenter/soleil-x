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

local F = {}
package.loaded["ebb.src.functions"] = F

local T                 = require 'ebb.src.types'
local Util              = require 'ebb.src.util'

local Pre               = require 'ebb.src.prelude'
local R                 = require 'ebb.src.relations'
local specialization    = require 'ebb.src.specialization'
local semant            = require 'ebb.src.semant'
local phase             = require 'ebb.src.phase'
local M                 = require "ebb.src.main"

local function shallowcopy_table(tbl)
  local x = {}
  for k,v in pairs(tbl) do x[k] = v end
  return x
end

-------------------------------------------------------------------------------

local Function    = {}
Function.__index  = Function
F.Function        = Function
local function is_function(obj) return getmetatable(obj) == Function end
F.is_function     = is_function


-------------------------------------------------------------------------------
--[[ UserFunc:                                                             ]]--
-------------------------------------------------------------------------------


function F.NewFunction(func_ast, luaenv)
  local special = specialization.specialize(luaenv, func_ast)

  local ufunc = setmetatable({
    _decl_ast     = special,
    _versions     = {}, -- the versions table is nested
    _name         = special.id,
  }, Function)

  M.decls():insert(M.AST.NewFunction(ufunc))

  return ufunc
end

function Function:Name()
  return self._name
end

function Function:setname(name)
  if type(name) ~= 'string' then error('expected string argument', 2) end
  if self._typed_at_least_once then
    error('cannot re-name a function after it has been compiled once', 2)
  end
  self._name = name
  self._decl_ast.id = name
end

local get_ufunc_typetable =
Util.memoize_from(2, function(calldepth, ufunc, relset, ...)
  calldepth = calldepth+1 -- account for the memoization wrapper
  -- ... are string arguments to function call

  local relation      = R.is_subset(relset) and relset:Relation() or relset
  -- make a safe copy that we can explicitly type annotate
  local aname_ast     = ufunc._decl_ast:alpha_rename()

  -- process the first argument's type annotation.  Consistent? Present?
  local annotation    = aname_ast.ptypes[1]
  if annotation then
    local arel = annotation.relation
    if arel ~= relation then
      error('The supplied relation did not match the parameter '..
            'annotation:\n  '..relation:Name()..' vs. '..arel:Name(),
            calldepth)
    end
  else
    -- add an annotation if none was present
    aname_ast.ptypes[1] = T.key(relation)
  end

  -- process the remaining arguments' type annotations.
  for i,str in ipairs({...}) do
    local annotation = aname_ast.ptypes[i+1]
    if annotation then
      error('Secondary string arguments to functions should be '..
            'untyped arguments', calldepth)
    end
    aname_ast.ptypes[i+1] = T.internal(str)
  end

  -- now actually type and phase check
  local typed_ast       = semant.check( aname_ast )
  local phase_results   = phase.phasePass( typed_ast )

  return {
    typed_ast       = typed_ast,
    field_use       = phase_results.field_use,
    global_use      = phase_results.global_use,
    inserts         = phase_results.inserts,
    deletes         = phase_results.deletes,
    versions        = terralib.newlist(),
  }
end)

function Function:_get_typechecked(calldepth, relset, strargs)
  return get_ufunc_typetable(calldepth+1, self, relset, unpack(strargs))
end

function Function:_Get_Type_Version_Table(calldepth, relset, ...)
  self._typed_at_least_once = true
  if not (R.is_subset(relset) or R.is_relation(relset)) then
    error('Functions must be executed over a relation or subset, but '..
          'argument was neither: '..tostring(relset), calldepth)
  end

  -- unpack direct arguments and/or launch parameters
  local args    = {...}
  local params  = {}
  if type(args[#args]) == 'table' then
    params = args[#args]
    args[#args] = nil
  end

  -- check that number of arguments matches, allowing for the
  -- extra first argument in the function signature that is a
  -- key for the relation being mapped over
  local narg_expected = #self._decl_ast.params - 1
  if narg_expected ~= #args then
    error('Function was expecting '..tostring(narg_expected)..
          ' arguments, but got '..tostring(#args), calldepth)
  end
  -- Also, right now we restrict all secondary arguments to be strings
  for i,a in ipairs(args) do
    if type(a) ~= 'string' then
      error('Argument '..tostring(i)..' was expected to be a string; '..
            'Secondary arguments to functions mapped over relations '..
            'must be strings.', calldepth)
    end
  end
  if self._decl_ast.exp then
    error('Functions executed over relations should not return values',
          calldepth)
  end

  -- get the appropriately typed version of the function
  -- and a collection of all the versions associated with it...
  local typeversion = self:_get_typechecked(calldepth+1, relset, args)

  return typeversion
end

local function get_ufunc_version(ufunc, typeversion_table, relset, params)
  params          = params or {}
  local proc      = params.location or Pre.default_processor
  local relation  = R.is_subset(relset) and relset:Relation() or relset

  return {
    ufunc           = ufunc,
    typtable        = typeversion_table,
    relation        = relation,
    subset          = R.is_subset(relset) and relset or nil,
    proc            = proc,
  }
end

local function get_func_call_params_from_args(...)
  local N = select('#',...)
  local last_arg = N > 0 and select(N,...) or nil
  if type(last_arg) == 'table' then return last_arg
                               else return {} end
end
function Function:_doForEach(relset, ...)
  local params      = get_func_call_params_from_args(...)
  local typeversion = self:_Get_Type_Version_Table(4, relset, ...)
  -- now we either retrieve or construct the appropriate function version
  local v = get_ufunc_version(self, typeversion, relset, params)
  M.stmts():insert(M.AST.ForEach(v.ufunc, v.relation, v.subset))
end
function Function:doForEach(relset, ...)
  self:_doForEach(relset, ...)
end
