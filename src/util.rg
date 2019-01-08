import 'regent'

local Exports = {}

local C = terralib.includecstring([[
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
]])

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

-- int, int -> T*
function Exports.range(first, last)
  local res = terralib.newlist()
  for i = first,last do
    res:insert(i)
  end
  return res
end

-------------------------------------------------------------------------------
-- Sets
-------------------------------------------------------------------------------

-- set(T) -> T*
function Exports.setToList(set)
  local list = terralib.newlist()
  for e,_ in pairs(set) do
    list:insert(e)
  end
  return list
end

-- T* -> set(T)
function Exports.listToSet(list)
  local set = {}
  for _,e in ipairs(list) do
    set[e] = true
  end
  return set
end

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

-- set(T), T* -> set(T)
function Exports.setSubList(set, list)
  local res = Exports.copyTable(set)
  for _,e in ipairs(list) do
    res[e] = nil
  end
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

-- terralib.struct -> string*
function Exports.fieldNames(s)
  return s.entries:map(function(e)
    local name,typ = Exports.parseStructEntry(e)
    return name
  end)
end

-- string, terralib.struct, string*, map(string,terralib.type)?
--   -> terralib.struct
function Exports.deriveStruct(newName, origStruct, keptFlds, addedFlds)
  addedFlds = addedFlds or {}
  local origFlds = {}
  for _,e in ipairs(origStruct.entries) do
    local name,typ = Exports.parseStructEntry(e)
    assert(not addedFlds[name])
    origFlds[name] = typ
  end
  local newStruct = terralib.types.newstruct(newName)
  for _,name in ipairs(keptFlds) do
    local typ = assert(origFlds[name])
    newStruct.entries:insert({name,typ})
  end
  for name,typ in pairs(addedFlds) do
    newStruct.entries:insert({name,typ})
  end
  return newStruct
end

-------------------------------------------------------------------------------
-- Filesystem
-------------------------------------------------------------------------------

terra Exports.createDir(name : &int8)
  var mode = 493 -- octal 0755 = rwxr-xr-x
  var res = C.mkdir(name, mode)
  if res < 0 then
    var stderr = C.fdopen(2, 'w')
    C.fprintf(stderr, 'Cannot create directory %s: ', name)
    C.fflush(stderr)
    C.perror('')
    C.fflush(stderr)
    C.exit(1)
  end
end

terra Exports.openFile(name : &int8, mode : &int8) : &C.FILE
  var file = C.fopen(name, mode)
  if file == nil then
    var stderr = C.fdopen(2, 'w')
    C.fprintf(stderr, 'Cannot open file %s in mode "%s": ', name, mode)
    C.fflush(stderr)
    C.perror('')
    C.fflush(stderr)
    C.exit(1)
  end
  return file
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
-- Regent metaprogramming
-------------------------------------------------------------------------------

-- regentlib.symbol, int, regentlib.rexpr, terralib.type -> regentlib.rquote
function Exports.emitRegionTagAttach(r, tag, value, typ)
  return rquote
    var info : typ[1]
    info[0] = value
    regentlib.c.legion_logical_region_attach_semantic_information(
      __runtime(), __raw(r), tag, [&typ](info), [sizeof(typ)], false)
  end
end

-- intXd, intXd, terralib.struct -> regentlib.task
function Exports.mkPartitionEqually(r_istype, cs_istype, fs)
  local partitionEqually
  if r_istype == int3d and cs_istype == int3d then
    __demand(__inline)
    task partitionEqually(r : region(ispace(int3d), fs),
                          cs : ispace(int3d),
                          halo : int3d,
                          offset : int3d)
      var Nx = r.bounds.hi.x - 2*halo.x + 1; var ntx = cs.bounds.hi.x + 1
      var Ny = r.bounds.hi.y - 2*halo.y + 1; var nty = cs.bounds.hi.y + 1
      var Nz = r.bounds.hi.z - 2*halo.z + 1; var ntz = cs.bounds.hi.z + 1
      regentlib.assert(Nx % ntx == 0, "Uneven partitioning on x")
      regentlib.assert(Ny % nty == 0, "Uneven partitioning on y")
      regentlib.assert(Nz % ntz == 0, "Uneven partitioning on z")
      regentlib.assert(-ntx <= offset.x and offset.x <= ntx, "offset.x too large")
      regentlib.assert(-nty <= offset.y and offset.y <= nty, "offset.y too large")
      regentlib.assert(-ntz <= offset.z and offset.z <= ntz, "offset.z too large")
      var coloring = regentlib.c.legion_domain_point_coloring_create()
      for c_real in cs do
        var c = (c_real - offset + {ntx,nty,ntz}) % {ntx,nty,ntz}
        var rect = rect3d{
          lo = int3d{halo.x + (Nx/ntx)*(c.x),
                     halo.y + (Ny/nty)*(c.y),
                     halo.z + (Nz/ntz)*(c.z)},
          hi = int3d{halo.x + (Nx/ntx)*(c.x+1) - 1,
                     halo.y + (Ny/nty)*(c.y+1) - 1,
                     halo.z + (Nz/ntz)*(c.z+1) - 1}}
        if c.x == 0 then rect.lo.x -= halo.x end
        if c.y == 0 then rect.lo.y -= halo.y end
        if c.z == 0 then rect.lo.z -= halo.z end
        if c.x == ntx-1 then rect.hi.x += halo.x end
        if c.y == nty-1 then rect.hi.y += halo.y end
        if c.z == ntz-1 then rect.hi.z += halo.z end
        regentlib.c.legion_domain_point_coloring_color_domain(coloring, c_real, rect)
      end
      var p = partition(disjoint, r, coloring, cs)
      regentlib.c.legion_domain_point_coloring_destroy(coloring)
      return p
    end
  elseif r_istype == int2d and cs_istype == int2d then
    __demand(__inline)
    task partitionEqually(r : region(ispace(int2d), fs),
                          cs : ispace(int2d),
                          halo : int2d,
                          offset : int2d)
      var Nx = r.bounds.hi.x - 2*halo.x + 1; var ntx = cs.bounds.hi.x + 1
      var Ny = r.bounds.hi.y - 2*halo.y + 1; var nty = cs.bounds.hi.y + 1
      regentlib.assert(Nx % ntx == 0, "Uneven partitioning on x")
      regentlib.assert(Ny % nty == 0, "Uneven partitioning on y")
      regentlib.assert(-ntx <= offset.x and offset.x <= ntx, "offset.x too large")
      regentlib.assert(-nty <= offset.y and offset.y <= nty, "offset.y too large")
      var coloring = regentlib.c.legion_domain_point_coloring_create()
      for c_real in cs do
        var c = (c_real - offset + {ntx,nty}) % {ntx,nty}
        var rect = rect2d{
          lo = int2d{halo.x + (Nx/ntx)*(c.x),
                     halo.y + (Ny/nty)*(c.y)},
          hi = int2d{halo.x + (Nx/ntx)*(c.x+1) - 1,
                     halo.y + (Ny/nty)*(c.y+1) - 1}}
        if c.x == 0 then rect.lo.x -= halo.x end
        if c.y == 0 then rect.lo.y -= halo.y end
        if c.x == ntx-1 then rect.hi.x += halo.x end
        if c.y == nty-1 then rect.hi.y += halo.y end
        regentlib.c.legion_domain_point_coloring_color_domain(coloring, c_real, rect)
      end
      var p = partition(disjoint, r, coloring, cs)
      regentlib.c.legion_domain_point_coloring_destroy(coloring)
      return p
    end
  elseif r_istype == int1d and cs_istype == int3d then
    __demand(__inline)
    task partitionEqually(r : region(ispace(int1d), fs),
                          cs : ispace(int3d),
                          halo : int64,
                          offset : int3d)
      var N = int64(r.bounds.hi - 2*halo + 1)
      var ntx = cs.bounds.hi.x + 1
      var nty = cs.bounds.hi.y + 1
      var ntz = cs.bounds.hi.z + 1
      regentlib.assert(N % (ntx*nty*ntz) == 0, "Uneven partitioning")
      regentlib.assert(-ntx <= offset.x and offset.x <= ntx, "offset.x too large")
      regentlib.assert(-nty <= offset.y and offset.y <= nty, "offset.y too large")
      regentlib.assert(-ntz <= offset.z and offset.z <= ntz, "offset.z too large")
      var coloring = regentlib.c.legion_domain_point_coloring_create()
      for c_real in cs do
        var c = (c_real - offset + {ntx,nty,ntz}) % {ntx,nty,ntz}
        var rect = rect1d{
          lo = halo + (N/ntx/nty/ntz)*(c.x*nty*ntz+c.y*ntz+c.z),
          hi = halo + (N/ntx/nty/ntz)*(c.x*nty*ntz+c.y*ntz+c.z+1) - 1}
        if c.x == 0 and
           c.y == 0 and
           c.z == 0 then rect.lo -= halo end
        if c.x == ntx-1 and
           c.y == nty-1 and
           c.z == ntz-1 then rect.hi += halo end
        regentlib.c.legion_domain_point_coloring_color_domain(coloring, c_real, rect)
      end
      var p = partition(disjoint, r, coloring, cs)
      regentlib.c.legion_domain_point_coloring_destroy(coloring)
      return p
    end
  else assert(false) end
  return partitionEqually
end

-------------------------------------------------------------------------------
-- Error handling
-------------------------------------------------------------------------------

-- regentlib.rexpr, string, regentlib.rexpr* -> regentlib.rquote
function Exports.emitAssert(cond, format, ...)
  local args = terralib.newlist{...}
  return rquote
    if not cond then
      var stderr = C.fdopen(2, 'w')
      C.fprintf(stderr, format, [args])
      C.fflush(stderr)
      C.exit(1)
    end
  end
end

-------------------------------------------------------------------------------

return Exports
