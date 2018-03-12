import 'regent'

local Exports = {}

local C = regentlib.c
local UNIX = terralib.includecstring([[
#include <sys/stat.h>
#include <sys/types.h>
]])

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

-- map(K,V) -> K*
function Exports.keys(tbl)
  local res = terralib.newlist()
  for k,v in pairs(tbl) do
    res:insert(k)
  end
  return res
end

-- T*, (T -> bool) -> bool
function Exports.all(list, pred)
  assert(terralib.israwlist(list))
  for _,x in ipairs(list) do
    if not pred(x) then return false end
  end
  return true
end

-- T*, (T -> bool) -> bool
function Exports.any(list, pred)
  assert(terralib.israwlist(list))
  for _,x in ipairs(list) do
    if pred(x) then return true end
  end
  return false
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

-- (T -> bool) -> bool
function TerraList:all(pred)
  return Exports.all(self, pred)
end

-- (T -> bool) -> bool
function TerraList:any(pred)
  return Exports.any(self, pred)
end

-- () -> T*
function TerraList:copy()
  return self:map(function(x) return x end)
end

-- T -> int?
function TerraList:find(x)
  for i,elem in ipairs(self) do
    if elem == x then
      return i
    end
  end
  return nil
end

-- () -> set(T)
function TerraList:toSet()
  local set = {}
  for _,elem in ipairs(self) do
    set[elem] = true
  end
  return set
end

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

-- () -> T*
function TerraList:reverse()
  local res = terralib.newlist()
  for i = #self, 1, -1 do
    res:insert(self[i])
  end
  return res
end

-- () -> T
function TerraList:pop()
  local res = self[#self]
  self[#self] = nil
  return res
end

-------------------------------------------------------------------------------
-- Strings
-------------------------------------------------------------------------------

-- string -> string*
function string:split(sep)
  local fields = terralib.newlist()
  local pattern = string.format("([^%s]+)", sep)
  self:gsub(pattern, function(c) fields:insert(c) end)
  return fields
end

-- string -> bool
function string:startswith(subStr)
  return self:sub(1, subStr:len()) == subStr
end

-- string -> bool
function string:endswith(subStr)
  return self:sub(self:len() - subStr:len() + 1, self:len()) == subStr
end

terra Exports.concretize(str : &int8) : int8[256]
  var res : int8[256]
  C.strncpy(&[&int8](res)[0], str, [uint64](256))
  [&int8](res)[255] = [int8](0)
  return res
end

-------------------------------------------------------------------------------
-- Terra type helpers
-------------------------------------------------------------------------------

-- {string,terralib.type} | {field:string,type:terralib.type} ->
--   string, terralib.type
function Exports.parseStructEntry(entry)
  if terralib.israwlist(entry) and #entry == 2 then
    return entry[1], entry[2]
  elseif entry.field and entry.type then
    return entry.field, entry.type
  else assert(false) end
end

-- terralib.type -> string, string
local function cTypeDecl(typ)
  if typ:isarray() then
    local decl, mods = cTypeDecl(typ.type)
    return decl, '['..tostring(typ.N)..']'..mods
  elseif typ == int    then return 'int',      ''
  elseif typ == int8   then return 'int8_t',   ''
  elseif typ == int16  then return 'int16_t',  ''
  elseif typ == int32  then return 'int32_t',  ''
  elseif typ == int64  then return 'int64_t',  ''
  elseif typ == uint   then return 'unsigned', ''
  elseif typ == uint8  then return 'uint8_t',  ''
  elseif typ == uint16 then return 'uint16_t', ''
  elseif typ == uint32 then return 'uint32_t', ''
  elseif typ == uint64 then return 'uint64_t', ''
  elseif typ == bool   then return 'bool',     ''
  elseif typ == float  then return 'float',    ''
  elseif typ == double then return 'double',   ''
  else assert(false) end
end

-- terralib.struct, bool?, string? -> string
function Exports.prettyPrintStruct(s, cStyle, indent)
  indent = indent or ''
  local lines = terralib.newlist()
  if s.name:startswith('anon') then
    lines:insert('struct {')
  else
    lines:insert('struct '..s.name..' {')
  end
  for _,e in ipairs(s.entries) do
    local name, typ = Exports.parseStructEntry(e)
    local s1 = cStyle and '' or (name..' : ')
    local s3 = cStyle and (' '..name) or ''
    local s2, s4 = '', ''
    if typ:isstruct() then
      s2 = Exports.prettyPrintStruct(typ, cStyle, indent..'  ')
    elseif cStyle then
      s2, s4 = cTypeDecl(typ)
    else
      s2 = tostring(typ)
    end
    lines:insert(indent..'  '..s1..s2..s3..s4..';')
  end
  lines:insert(indent..'}')
  return lines:join('\n')
end

-------------------------------------------------------------------------------
-- Filesystem
-------------------------------------------------------------------------------

terra Exports.mkdir(name : rawstring)
  var mode = 493 -- octal 0755 = rwxr-xr-x
  var res = UNIX.mkdir(name, mode);
  if res < 0 then
    var stderr = C.fdopen(2, 'w')
    C.fprintf(stderr, 'Cannot create directory %s: ', name)
    C.fflush(stderr)
    C.perror('')
    C.fflush(stderr)
    C.exit(1)
  end
end

-------------------------------------------------------------------------------

return Exports
