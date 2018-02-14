local Exports = {}

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

-------------------------------------------------------------------------------
-- Structs
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

-------------------------------------------------------------------------------

return Exports
