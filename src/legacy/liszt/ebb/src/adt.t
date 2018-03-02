--[[
  ADT - Algebraic Data Types
    Language for Abstract Syntax Trees and other data described by
    Algebraic Data Types in DSL development

    Based off of Zephyr from Python.

    Original Lua translation by Zach Devito (ASDL in Terra)
    Modified and Customized by Gilbert Bernstein (renamed ADT)
    March 2016
--]]
local ADT = require('ebb.src.adt_impl')

local function repack_ctxt(ctxt)
  local tbl = {}
  for k,v in pairs(ctxt.definitions) do tbl[k] = v end
  return tbl
end

local function commonExprStmt(self, lexer, namespace_name)
  local defs, externs = ADT._INTERNAL_terra_dsl_parse(lexer)

  local function constructor(env_fn)
    local ctxt = ADT.newcontext()
    ctxt:TerraDSLDefine(defs, externs, env_fn())
    return repack_ctxt(ctxt)
  end

  return constructor
end
local function handleExpression(self, lexer)
  lexer:expect('ADT')
  return commonExprStmt(self, lexer)
end
local function handleStatement(self, lexer)
  lexer:expect('ADT')
  local namespace_name  = lexer:expect(lexer.name).value
  local assigntuple     = { namespace_name }

  return commonExprStmt(self, lexer, namespace_name), assigntuple
end

local adt_language = {
  name          = 'adt_language',
  entrypoints   = {'ADT'},
  keywords      = {
    'attributes', 'extern'
  },

  expression      = handleExpression,
  statement       = handleStatement,
  localstatement  = handleStatement,
}

return adt_language