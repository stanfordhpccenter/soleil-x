--[[
  ADT - Algebraic Data Types
    Language for Abstract Syntax Trees and other data described by
    Algebraic Data Types in DSL development

    Based off of Zephyr from Python.

    Original Lua translation by Zach Devito (ASDL in Terra)
    Modified and Customized by Gilbert Bernstein (renamed ADT)
    March 2016
--]]

local newlist, islist = terralib.newlist, terralib.islist

--[[
form of parse tree:
ROOT            = { <def>* }
<def>           = { name      = <Ident>
                    type      = <type> }
<type>          = <product>
                | <sum>
<sum>           = { kind      = "sum"
                    constructors = { <constructor>* } }
<constructor>   = { name      = <Ident>
                    fields    = { <field>* } }
<product>       = { kind      = "product"
                    fields    = { <field>* } }
<field>         = { type      = <Ident>
                    optional  = <bool>
                    list      = <bool>
                    name      = <Ident> }
<Ident>         = [a-zA-Z_][a-zA-Z0-9_]*
--]]

local function parse_from_terra(lexer)
  local function parseField()
    local f = {}
    f.name  = lexer:expect(lexer.name).value
              lexer:expect(':')
    f.type  = lexer:expect(lexer.name).value
    if        lexer:nextif("?") then f.optional  = true
    elseif    lexer:nextif("*") then f.list      = true   end
    return f
  end
  local function parseFields()
    local fields = newlist()
    lexer:expect('{')
    if not lexer:matches('}') then repeat
      fields:insert( parseField() )
    until not lexer:nextif(',') end
    lexer:expect('}')
    return fields
  end
  local function parseProduct()
    local p   = {}
    p.kind    = "product"
    p.fields  = parseFields()
    return p
  end

  local function parseConstructor()
    local c   = {}
    c.name    = lexer:expect(lexer.name).value
    if lexer:matches('{') then
      c.fields  = parseFields()
    end
    return c
  end
  local function parseSum()
    local sum         = {}
    sum.kind          = "sum"
    sum.constructors  = newlist()
    repeat
      sum.constructors:insert( parseConstructor() )
    until not lexer:nextif("|")
    if lexer:nextif("attributes") then
      local attributes = parseFields()
      for i,ctor in ipairs(sum.constructors) do
        ctor.fields = ctor.fields or newlist()
        for i,a in ipairs(attributes) do
          ctor.fields:insert(a)
        end
      end
    end
    return sum
  end

  local function parseType()
    if lexer:matches('{') then return parseProduct()
                          else return parseSum()      end
  end
  local function parseExtern()
    local x = {}
    x.name  = lexer:expect(lexer.name).value
    if lexer:matches('function') then
      x.funcexpr  = lexer:luaexpr()
    else
      x.funcname  = lexer:expect(lexer.name).value
      lexer:ref(x.funcname)
    end
    return x
  end
  local function parseDefinitions()
    local ds = newlist()
    local xs = newlist()
    while not lexer:matches('end') do
      if lexer:nextif('extern') then
        xs:insert(parseExtern())
      else
        local d = {}
        d.name  = lexer:expect(lexer.name).value
                  lexer:expect('=')
        d.type  = parseType()
        ds:insert(d)
      end
    end
    return ds, xs
  end

  local defs, externs = parseDefinitions()
  lexer:expect("end")
  return defs, externs
end


local function build_check_builtin(t)
  return function(v) return type(v) == t end
end

local function build_check_optional(checkt)
  return function(v) return v == nil or checkt(v) end
end

local function build_check_list(checkt)
  return function(vs)
    if not islist(vs) then return false end
    for i,e in ipairs(vs) do
      if not checkt(e) then return false, i end
    end
    return true
  end
end

local defaultchecks = {}
for _,str in ipairs({
  'nil', 'number', 'string', 'boolean', 'table',
  'thread', 'userdata', 'cdata', 'function',
}) do
  defaultchecks[str] = build_check_builtin(str)
end
defaultchecks['any'] = function() return true end

local Context = {}
function Context:__index(idx)
    local d = self.definitions[idx]
    if d ~= nil then return d end
    return getmetatable(self)[idx]
end

local function newcontext()
  return setmetatable({
    checks        = setmetatable({},{__index = defaultchecks}),
    list          = {}, -- hold cached field checks
    optional      = {}, -- hold cached field checks
    definitions   = {},
  }, Context)
end

function Context:GetCheckForField(field)
  local ctxt = self
  -- gets the primitive or already-defined check function
  local check = ctxt.checks[field.type]
  if not check then error("type not defined: "..tostring(field.type), 5) end

  local function get(tbl,ctor)
    if not tbl[field.type] then
      tbl[field.type] = ctor(check)
    end
    return tbl[field.type]
  end
  if field.list then
    return get(ctxt.list,build_check_list)
  elseif field.optional then
    return get(ctxt.optional,build_check_optional)
  end
  return check
end

local function reporterr(i,typname,name,v,ii)
  local fmt = "bad argument #%d to '%s' expected '%s' but found '%s'"
  if ii then v,fmt = v[ii],fmt .. " at list index %d" end
  local err = string.format(fmt,i,name,typname,type(v),ii)
  local mt  = getmetatable(v)
  if mt then err = string.format("%s (metatable = %s)",err,tostring(mt)) end
  error(err, 3)
end
function Context:Extern(name,istype)
  self.checks[name] = istype
end

function Context:Define(text)
  local defs = parseAll(text)
  return self:Shared_Define_Code(defs)
end

function Context:TerraDSLDefine(defs, externs, env)
  for _,x in ipairs(externs) do
    local status, test_func
    if x.funcname then
      test_func = env[x.funcname]
      if not test_func then
        error('Could not find extern test function '..x.funcname..
              ' in the enclosing lua scope')
      end
    elseif x.funcexpr then
      status, test_func = pcall(function()
        return x.funcexpr(env)
      end)
      if not status then
        error("error while evaluating lua expression for "..
              "extern '"..x.name.."'")
      end
    end
    if type(test_func) ~= 'function' then
      error("expected 'extern "..x.name.."' to specify a function")
    end
    self:Extern(x.name,test_func)
  end
  return self:Shared_Define_Code(defs)
end

function Context:Shared_Define_Code(defs)
  local ctxt = self

  -- Declaration of Classes
  --[[
    Called for every Sum, Product, and Case Class
    Post-Condition:
      1)  definitions contains an entry { members = {} } for each class
          this SET will be filled out later so that...
      2)  a CHECK function is defined for each class that checks whether the
          argument is a member class (i.e. exactly this class,
          or for a sum class, one of its case classes)
  --]]
  local function DeclareClass(name)
    assert(not ctxt.definitions[name],
           "class name "..name.." already defined")
    local m = {}
    ctxt.definitions[name] = { members = m }
    ctxt.checks[name] = function(v)
      return m[getmetatable(v) or false] or false
    end
  end
  for i,d in ipairs(defs) do
    DeclareClass(d.name)
    if d.type.kind == "sum" then
      for i,c in ipairs(d.type.constructors) do
        DeclareClass(c.name)
      end
    end
  end

  -- Definition of Classes
  local function DefineClass(name, fields)
    local mt    = {}
    local class = ctxt.definitions[name]

    if fields then
      local names     = newlist()
      local checks    = newlist()
      local typnames  = newlist()
      for i,f in ipairs(fields) do
        names:insert(f.name)
        typnames:insert(f.list and f.type.."*" or f.type)
        checks:insert(ctxt:GetCheckForField(f))
      end

      -- Constructor
      function mt:__call(...)
        local obj = {}
        for i=1, #names do
          local v = select(i,...)
          local c,ii = checks[i](v)
          if not c then reporterr(i,typnames[i],name,v,ii) end
          obj[names[i]] = v
        end
        return setmetatable(obj,self)
      end
      -- Print
      function class:__tostring()
        local members = newlist()
        for i,f in ipairs(fields) do
          local v = self[f.name]
          if f.optional and v == nil then -- do nothing
          else
            if f.list then
              local elems = newlist()
              for i,e in ipairs(v) do elems:insert(tostring(e)) end
              members:insert(f.name.." = {"..elems:concat(",").."}")
            else
              members:insert(f.name.." = "..tostring(v))
            end
          end
        end
        return name..'('..members:concat(",")..')'
      end
    else
      function class:__tostring() return name end
    end
    function mt:__tostring() return "Class("..name..")" end
    function mt:__newindex(k,v)
      for c,_ in pairs(self.members) do rawset(c,k,v) end
    end
    local check = assert(ctxt.checks[name])
    function class.check(obj)
      return check(obj)
    end
    class.__index = class
    class.members[class] = true
    setmetatable(class,mt)
    return class
  end
  for i,d in ipairs(defs) do
    if d.type.kind == "sum" then
      local parent = DefineClass(d.name, nil)
      for i,c in ipairs(d.type.constructors) do
        local child = DefineClass(c.name, c.fields)
        parent.members[child] = true -- mark subclass as member of parent
        child.kind = c.name
        if not c.fields then -- 0-ary constructor value
          ctxt.definitions[c.name] = setmetatable({}, child)
        end
      end
    else
      DefineClass(d.name, d.type.fields)
    end
  end
end

return { newcontext = newcontext,
         _INTERNAL_terra_dsl_parse = parse_from_terra }