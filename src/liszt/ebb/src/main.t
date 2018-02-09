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

import 'ebb.src.adt'
import 'regent'

local M = {}
package.loaded['ebb.src.main'] = M

local L = require 'ebblib'
local isField    = function(x) return L.is_field(x) end
local isFunction = function(x) return L.is_function(x) end
local isGlobal   = function(x) return L.is_global(x) end
local isRelation = function(x) return L.is_relation(x) end
local isSubset   = function(x) return L.is_subset(x) end

local RG  = regentlib

-------------------------------------------------------------------------------
-- Control language AST
-------------------------------------------------------------------------------

-- ExprConst = boolean | number | ExprConst[N]

-- any -> boolean
local function isExprConst(x)
  if type(x) == 'boolean' or type(x) == 'number' then
    return true
  end
  if not terralib.israwlist(x) then
    return false
  end
  for _,elem in ipairs(x) do
    if not isExprConst(elem) then
      return false
    end
  end
  return true
end
M.isExprConst = isExprConst

local ADT AST
  Decl = NewField { fld : Field }
       | NewFunction { fun : Function }
       | NewGlobal { global : Global, init : Expr }
       | NewRelation { rel : Relation }
  Stmt = Block { stmts : Stmt* }
       | ForEach { fun : Function, rel : Relation, subset : Subset? }
       | If { cond : Cond, thenBlock : Stmt?, elseBlock : Stmt? }
       | FillField { fld : Field, val : Expr }
       | SetGlobal { global : Global, expr : Expr }
       | While { cond : Cond, spmd : boolean, body : Stmt? }
       | Do { spmd : boolean, body : Stmt? }
       | Print { fmt : string, globals : Global* }
       | Dump { rel : Relation, flds : string*, file : string, vals : Expr* }
       | Load { rel : Relation, flds : string*, file : string, vals : Expr* }
       | Inline { quot : RQuote }
       | Error { msg : string }
  Cond = Literal { val : boolean }
       | And { lhs : Cond, rhs : Cond }
       | Or { lhs : Cond, rhs : Cond }
       | Not { cond : Cond }
       | Compare { op : string, lhs : Expr, rhs : Expr }
  Expr = Const { val : ExprConst }
       | GetGlobal { global : Global }
       | BinaryOp { op : string, lhs : Expr, rhs : Expr }
       | UnaryOp { op : string, arg : Expr }
       | Array { elems : Expr* }
       | Index { base : Expr, index : number }
       | ReadConfig { name : string }
       | Cond2Expr { cond : Cond }
  extern ExprConst isExprConst
  extern Field     isField
  extern Function  isFunction
  extern Global    isGlobal
  extern Relation  isRelation
  extern Subset    isSubset
  extern RQuote    function(x) return RG.is_rquote(x) end
end
M.AST = AST

-------------------------------------------------------------------------------
-- API call logs
-------------------------------------------------------------------------------

local decls = terralib.newlist() -- AST.Decl*

-- () -> AST.Decl*
function M.decls()
  return decls
end

local scopes = terralib.newlist({terralib.newlist()}) -- (AST.Stmt*)*
local stack = terralib.newlist() -- (AST.If | AST.While)*

-- () -> AST.Stmt*
function M.stmts()
  return scopes[#scopes]
end

-------------------------------------------------------------------------------
-- Pretty-printing AST nodes
-------------------------------------------------------------------------------

function AST.Expr:__tostring()
  error('Abstract Method')
end
function AST.Const:__tostring()
  local str
  if type(self.val) == 'table' then
    local elemStrs =
      terralib.newlist(self.val):map(function(x) return tostring(x) end)
    str = '['..elemStrs:join(',')..']'
  else
    str = tostring(self.val)
  end
  return 'Const('..str..')'
end
function AST.GetGlobal:__tostring()
  return 'GetGlobal('..self.global:Name()..')'
end
function AST.BinaryOp:__tostring()
  return
    'BinaryOp('..self.op..','..tostring(self.lhs)..','..tostring(self.rhs)..')'
end
function AST.UnaryOp:__tostring()
  return 'UnaryOp('..self.op..','..tostring(self.arg)..')'
end
function AST.Array:__tostring()
  local elemStrs = self.elems:map(function(x) return tostring(x) end)
  return 'Array('..elemStrs:join(',')..')'
end
function AST.Index:__tostring()
  return 'Index('..tostring(self.base)..','..tostring(self.index)..')'
end
function AST.ReadConfig:__tostring()
  return 'ReadConfig('..self.name..')'
end
function AST.Cond2Expr:__tostring()
  return 'Cond2Expr(...)'
end

-------------------------------------------------------------------------------
-- Operations on AST nodes
-------------------------------------------------------------------------------

-- string, ExprConst | AST.Expr, ExprConst | AST.Expr -> AST.Expr
local function implBinaryOp(op, lhs, rhs)
  if isExprConst(lhs) then lhs = AST.Const(lhs) end
  if isExprConst(rhs) then rhs = AST.Const(rhs) end
  return AST.BinaryOp(op, lhs, rhs)
end

-- ExprConst | AST.Expr, ExprConst | AST.Expr -> AST.Expr
function AST.Expr.__add(lhs, rhs) return implBinaryOp('+', lhs, rhs) end
function AST.Expr.__sub(lhs, rhs) return implBinaryOp('-', lhs, rhs) end
function AST.Expr.__mul(lhs, rhs) return implBinaryOp('*', lhs, rhs) end
function AST.Expr.__div(lhs, rhs) return implBinaryOp('/', lhs, rhs) end
function AST.Expr.__mod(lhs, rhs) return implBinaryOp('%', lhs, rhs) end

-- AST.Expr -> AST.Expr
function AST.Expr.__unm(x)
  return AST.UnaryOp('-', x)
end

local function wrapIndex(class)
  function class.__index(obj, key)
    if type(key) == 'number' and key == math.floor(key) then
      return AST.Index(obj, key)
    end
    local val = rawget(class, key)
    if val then return val end
    error('Key '..key..' not found in M.AST.Expr subclass')
  end
end
wrapIndex(AST.Const)
wrapIndex(AST.GetGlobal)
wrapIndex(AST.BinaryOp)
wrapIndex(AST.UnaryOp)
wrapIndex(AST.Array)
wrapIndex(AST.Index)
wrapIndex(AST.ReadConfig)

-- ExprConst | AST.Expr, ExprConst | AST.Expr -> AST.Expr
function M.MAX(lhs, rhs)
  if isExprConst(lhs) then lhs = AST.Const(lhs) end
  if isExprConst(rhs) then rhs = AST.Const(rhs) end
  return AST.BinaryOp('max', lhs, rhs)
end

-- ExprConst | AST.Expr, ExprConst | AST.Expr -> AST.Expr
function M.MIN(lhs, rhs)
  if isExprConst(lhs) then lhs = AST.Const(lhs) end
  if isExprConst(rhs) then rhs = AST.Const(rhs) end
  return AST.BinaryOp('min', lhs, rhs)
end

-- boolean | AST.Cond, boolean | AST.Cond -> AST.Cond
function M.AND(lhs, rhs)
  if type(lhs) == 'boolean' then lhs = AST.Literal(lhs) end
  if type(rhs) == 'boolean' then rhs = AST.Literal(rhs) end
  return AST.And(lhs, rhs)
end

-- boolean | AST.Cond, boolean | AST.Cond -> AST.Cond
function M.OR(lhs, rhs)
  if type(lhs) == 'boolean' then lhs = AST.Literal(lhs) end
  if type(rhs) == 'boolean' then rhs = AST.Literal(rhs) end
  return AST.Or(lhs, rhs)
end

-- boolean | AST.Cond -> AST.Cond
function M.NOT(cond)
  if type(cond) == 'boolean' then cond = AST.Literal(cond) end
  return AST.Not(cond)
end

-- string, ExprConst | AST.Expr, ExprConst | AST.Expr -> AST.Cond
local function implCompare(op, lhs, rhs)
  if isExprConst(lhs) then lhs = AST.Const(lhs) end
  if isExprConst(rhs) then rhs = AST.Const(rhs) end
  return AST.Compare(op, lhs, rhs)
end

-- ExprConst | AST.Expr, ExprConst | AST.Expr -> AST.Cond
function M.EQ(lhs, rhs) return implCompare('==', lhs, rhs) end
function M.NE(lhs, rhs) return implCompare('~=', lhs, rhs) end
function M.GT(lhs, rhs) return implCompare('>',  lhs, rhs) end
function M.LT(lhs, rhs) return implCompare('<',  lhs, rhs) end
function M.GE(lhs, rhs) return implCompare('>=', lhs, rhs) end
function M.LE(lhs, rhs) return implCompare('<=', lhs, rhs) end

-- (ExprConst | AST.Expr)* -> AST.Expr
function M.ARRAY(elems)
  return AST.Array(terralib.newlist(elems):map(function(e)
    return isExprConst(e) and AST.Const(e) or e
  end))
end

-- AST.Cond -> AST.Expr
function M.COND2EXPR(cond)
  return AST.Cond2Expr(cond)
end

-------------------------------------------------------------------------------
-- Control flow
-------------------------------------------------------------------------------

-- boolean | AST.Cond -> ()
function M.IF(cond)
  if type(cond) == 'boolean' then cond = AST.Literal(cond) end
  stack:insert(AST.If(cond, nil, nil))
  scopes:insert(terralib.newlist())
end

-- () -> ()
function M.ELSE()
  assert(#stack > 0)
  local wrapper = stack[#stack]
  assert(AST.If.check(wrapper))
  assert(not wrapper.thenBlock)
  wrapper.thenBlock = AST.Block(scopes[#scopes])
  scopes[#scopes] = terralib.newlist()
end

-- boolean | AST.Cond, boolean? -> ()
function M.WHILE(cond, spmd)
  spmd = spmd or false
  if type(cond) == 'boolean' then cond = AST.Literal(cond) end
  stack:insert(AST.While(cond, spmd, nil))
  scopes:insert(terralib.newlist())
end

-- boolean? -> ()
function M.DO(spmd)
  spmd = spmd or false
  stack:insert(AST.Do(spmd, nil))
  scopes:insert(terralib.newlist())
end

-- () -> ()
function M.END()
  assert(#stack > 0)
  local wrapper = stack[#stack]
  local stmts = scopes[#scopes]
  if AST.If.check(wrapper) then
    if wrapper.thenBlock then
      assert(not wrapper.elseBlock)
      wrapper.elseBlock = AST.Block(stmts)
    else
      wrapper.thenBlock = AST.Block(stmts)
    end
  elseif AST.While.check(wrapper) then
    wrapper.body = AST.Block(stmts)
  elseif AST.Do.check(wrapper) then
    wrapper.body = AST.Block(stmts)
  else assert(false) end
  scopes[#scopes-1]:insert(wrapper)
  scopes[#scopes] = nil
  stack[#stack] = nil
end

-------------------------------------------------------------------------------
-- Misc commands
-------------------------------------------------------------------------------

-- string, PRE.Global* -> ()
function M.PRINT(fmt, ...)
  M.stmts():insert(AST.Print(fmt, terralib.newlist({...})))
end

-- RG.rquote -> ()
function M.INLINE(quot)
  M.stmts():insert(AST.Inline(quot))
end

-- string -> ()
function M.ERROR(msg)
  M.stmts():insert(AST.Error(msg))
end
