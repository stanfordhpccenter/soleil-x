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


--[[
--
--  prelude.t
--
--    This file encapsulates assorted objects/concepts that will get
--    included into the public API.  Rather than pollute the ebblib.t
--    file specifying that interface in total, the implementations are
--    hidden here.
--
--]]

local Pre = {}
package.loaded["ebb.src.prelude"] = Pre
local M = require "ebb.src.main"

-------------------------------------------------------------------------------

local Global      = {}
Global.__index    = Global
Pre.Global        = Global
local function is_global(obj) return getmetatable(obj) == Global end
Pre.is_global     = is_global

local Constant    = {}
Constant.__index  = Constant
Pre.Constant      = Constant
local function is_constant(obj) return getmetatable(obj) == Constant end
Pre.is_constant   = is_constant

local Macro       = {}
Macro.__index     = Macro
Pre.Macro         = Macro
local function is_macro(obj) return getmetatable(obj) == Macro end
Pre.is_macro      = is_macro


-------------------------------------------------------------------------------

local T   = require 'ebb.src.types'

-------------------------------------------------------------------------------
--[[ Globals:                                                              ]]--
-------------------------------------------------------------------------------

function Pre.NewGlobal (name, typ, init)
  if not T.istype(typ) or not typ:isvalue() then
    error("Second argument to Global() must be an Ebb value type", 2)
  end
  if M.isExprConst(init) then
    -- Lift lua constant to AST node.
    if not T.luaValConformsToType(init, typ) then
      error("Third argument to Global() must be an "..
            "instance of type " .. tostring(typ), 2)
    end
    init = M.AST.Const(init)
  end
  local s  = setmetatable({_name=name, _type=typ}, Global)
  M.decls():insert(M.AST.NewGlobal(s, init))
  return s
end

function Global:__newindex(fieldname,value)
  error("Cannot assign members to Global object", 2)
end

function Global:set(val)
  if M.isExprConst(val) then
    -- Lift lua constant to AST node.
    if not T.luaValConformsToType(val, self._type) then
      error("value does not conform to type of global: "..
              tostring(self._type), 2)
    end
    val = M.AST.Const(val)
  end
  M.stmts():insert(M.AST.SetGlobal(self, val))
end

function Global:get()
  return M.AST.GetGlobal(self)
end

function Global:Name()
  return self._name
end

function Global:Type()
  return self._type
end

-------------------------------------------------------------------------------
--[[ Constants:                                                            ]]--
-------------------------------------------------------------------------------

local function deep_copy(tbl)
    if type(tbl) ~= 'table' then return tbl
    else
        local cpy = {}
        for i=1,#tbl do cpy[i] = deep_copy(tbl[i]) end
        return cpy
    end
end

function Pre.NewConstant (typ, init)
    if not T.istype(typ) or not typ:isvalue() then
        error("First argument to Constant() must be an "..
              "Ebb value type", 2)
    end
    if not T.luaValConformsToType(init, typ) then
        error("Second argument to Constant() must be a "..
              "value of type " .. tostring(typ), 2)
    end


    local c = setmetatable({_type=typ, _value=deep_copy(init)}, Constant)
    return c
end

function Constant:__newindex(fieldname,value)
  error("Cannot assign members to Constant object", 2)
end

function Constant:get()
  return deep_copy(self._value)
end

function Constant:Type()
  return deep_copy(self._value)
end

-------------------------------------------------------------------------------
--[[ LMacros:                                                              ]]--
-------------------------------------------------------------------------------

function Pre.NewMacro(generator)
    return setmetatable({genfunc=generator}, Macro)
end

function Macro:__newindex(fieldname,value)
  error("Cannot assign members to Macro object", 2)
end
