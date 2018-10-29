import 'regent'

local Exports = {}

-------------------------------------------------------------------------------
-- Numeric
-------------------------------------------------------------------------------

-- A -> bool
function Exports.isPosInt(x)
  return type(x) == 'number' and x == math.floor(x) and x > 0
end

-------------------------------------------------------------------------------
-- Tables
-------------------------------------------------------------------------------

-- table -> bool
function Exports.isEmpty(tbl)
  if not tbl then
    return true
  end
  for _,_ in pairs(tbl) do
    return false
  end
  return true
end

-- table -> table
function Exports.copyTable(tbl)
  local cpy = {}
  for k,v in pairs(tbl) do cpy[k] = v end
  return cpy
end

-- table -> int
function Exports.tableSize(tbl)
  local size = 0
  for _,_ in pairs(tbl) do
    size = size + 1
  end
  return size
end

-------------------------------------------------------------------------------
-- Lists
-------------------------------------------------------------------------------

local TerraList = getmetatable(terralib.newlist())

-- string? -> string
function TerraList:join(sep)
  sep = sep or ' '
  local res = ''
  for i,elem in ipairs(self) do
    if i > 1 then
      res = res..sep
    end
    res = res..tostring(elem)
  end
  return res
end

-- int, (() -> T) -> T*
function Exports.generate(n, generator)
  local res = terralib.newlist()
  for i = 1,n do
    res:insert(generator())
  end
  return res
end

-- () -> T*
function TerraList:flatten(res)
  res = res or terralib.newlist()
  for _,e in ipairs(self) do
    if terralib.israwlist(e) then
      e:flatten(res)
    else
      res:insert(e)
    end
  end
  return res
end

-------------------------------------------------------------------------------
-- Sets
-------------------------------------------------------------------------------

-- set(T) -> T
function Exports.setPop(set)
  local elem
  for e,_ in pairs(set) do
    elem = e
    break
  end
  set[elem] = nil
  return elem
end

-------------------------------------------------------------------------------
-- Strings
-------------------------------------------------------------------------------

-- string -> bool
function string:startswith(subStr)
  return self:sub(1, subStr:len()) == subStr
end

-- string -> bool
function string:endswith(subStr)
  return self:sub(self:len() - subStr:len() + 1, self:len()) == subStr
end

-------------------------------------------------------------------------------
-- Terra type helpers
-------------------------------------------------------------------------------

-- {string,terralib.type} | {field:string,type:terralib.type} ->
--   string, terralib.type
function Exports.parseStructEntry(entry)
  if terralib.israwlist(entry)
  and #entry == 2
  and type(entry[1]) == 'string'
  and terralib.types.istype(entry[2]) then
    return entry[1], entry[2]
  elseif type(entry) == 'table'
  and entry.field
  and entry.type then
    return entry.field, entry.type
  else assert(false) end
end

-- map(terralib.type,string)
local cBaseType = {
  [int]    = 'int',
  [int8]   = 'int8_t',
  [int16]  = 'int16_t',
  [int32]  = 'int32_t',
  [int64]  = 'int64_t',
  [uint]   = 'unsigned',
  [uint8]  = 'uint8_t',
  [uint16] = 'uint16_t',
  [uint32] = 'uint32_t',
  [uint64] = 'uint64_t',
  [bool]   = 'bool',
  [float]  = 'float',
  [double] = 'double',
}

-- terralib.type, bool, string -> string, string
local function typeDecl(typ, cStyle, indent)
  if typ:isarray() then
    local decl, mods = typeDecl(typ.type, cStyle, indent)
    decl = cStyle and decl or decl..'['..tostring(typ.N)..']'
    mods = cStyle and '['..tostring(typ.N)..']'..mods or mods
    return decl, mods
  elseif typ:isstruct() then
    if typ.name:startswith('anon') then
      return Exports.prettyPrintStruct(typ, cStyle, indent), ''
    elseif cStyle then
      return ('struct '..typ.name), ''
    else
      return typ.name, ''
    end
  elseif typ:isprimitive() then
    return (cStyle and cBaseType[typ] or tostring(typ)), ''
  else assert(false) end
end

-- terralib.struct, bool?, string? -> string
function Exports.prettyPrintStruct(s, cStyle, indent)
  indent = indent or ''
  local lines = terralib.newlist()
  local isUnion = false
  local entries = s.entries
  if #entries == 1 and terralib.israwlist(entries[1]) then
    isUnion = true
    entries = entries[1]
  end
  local name = s.name:startswith('anon') and '' or s.name
  local open =
    (cStyle and isUnion)         and ('union '..name..' {')          or
    (cStyle and not isUnion)     and ('struct '..name..' {')         or
    (not cStyle and isUnion)     and ('struct '..name..' { union {') or
    (not cStyle and not isUnion) and ('struct '..name..' {')         or
    assert(false)
  lines:insert(open)
  for _,e in ipairs(entries) do
    local name, typ = Exports.parseStructEntry(e)
    local s1 = cStyle and '' or (name..' : ')
    local s3 = cStyle and (' '..name) or ''
    local s2, s4 = typeDecl(typ, cStyle, indent..'  ')
    lines:insert(indent..'  '..s1..s2..s3..s4..';')
  end
  local close =
    (cStyle and isUnion)         and '}'   or
    (cStyle and not isUnion)     and '}'   or
    (not cStyle and isUnion)     and '} }' or
    (not cStyle and not isUnion) and '}'   or
    assert(false)
  lines:insert(indent..close)
  return lines:join('\n')
end

-------------------------------------------------------------------------------
-- Graphs
-------------------------------------------------------------------------------

-- Graph(T) = map(T,T*)

-- Graph(T) -> set(T)
function Exports.getNodes(graph)
  local nodes = {} -- set(T)
  for src,tgts in pairs(graph) do
    nodes[src] = true
    for _,tgt in ipairs(tgts) do
      nodes[tgt] = true
    end
  end
  return nodes
end

-- Graph(T) -> T*
function Exports.revTopoSort(graph)
  local unmarked = Exports.getNodes(graph) -- set(T)
  local tempMarked = {} -- set(T)
  local permMarked = {} -- set(T)
  local res = terralib.newlist() -- T*
  local function visit(src)
    if permMarked[src] then
      return true
    end
    if tempMarked[src] then
      return false
    end
    tempMarked[src] = true
    for _,tgt in ipairs(graph[src] or {}) do
      visit(tgt)
    end
    permMarked[src] = true
    res:insert(src)
    return true
  end
  while not Exports.isEmpty(unmarked) do
    if not visit(Exports.setPop(unmarked)) then
      return nil
    end
  end
  return res
end

-------------------------------------------------------------------------------

return Exports
