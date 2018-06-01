local C = terralib.includecstring [[
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
]]
local JSON = terralib.includec('json.h')
local UTIL = require 'util'

-------------------------------------------------------------------------------

-- NOTE: Type constructors are placed in the global namespace, so the schema
-- file can access them without imports.

local isStruct
local isSchemaT

-------------------------------------------------------------------------------

local StringMT = {}
StringMT.__index = StringMT

-- int -> String
function String(maxLen)
  assert(UTIL.isPosInt(maxLen))
  return setmetatable({
    maxLen = maxLen,
  }, StringMT)
end

-- SchemaT -> bool
local function isString(typ)
  return type(typ) == 'table' and getmetatable(typ) == StringMT
end

-------------------------------------------------------------------------------

local EnumMT = {}
EnumMT.__index = EnumMT

-- string* -> Enum
function Enum(...)
  local enum = {}
  for i,choice in ipairs({...}) do
    assert(type(choice) == 'string')
    enum[choice] = i-1
  end
  return setmetatable(enum, EnumMT)
end

-- SchemaT -> bool
local function isEnum(typ)
  return type(typ) == 'table' and getmetatable(typ) == EnumMT
end

-------------------------------------------------------------------------------

local ArrayMT = {}
ArrayMT.__index = ArrayMT

-- int, SchemaT -> Array
function Array(num, elemType)
  assert(UTIL.isPosInt(num))
  assert(isSchemaT(elemType))
  return setmetatable({
    num = num,
    elemType = elemType,
  }, ArrayMT)
end

-- SchemaT -> bool
local function isArray(typ)
  return type(typ) == 'table' and getmetatable(typ) == ArrayMT
end

-------------------------------------------------------------------------------

local UpToMT = {}
UpToMT.__index = UpToMT

-- int, SchemaT -> UpTo
function UpTo(max, elemType)
  assert(UTIL.isPosInt(max))
  assert(isSchemaT(elemType))
  return setmetatable({
    max = max,
    elemType = elemType,
  }, UpToMT)
end

-- SchemaT -> bool
local function isUpTo(typ)
  return type(typ) == 'table' and getmetatable(typ) == UpToMT
end

-------------------------------------------------------------------------------

local UnionMT = {}
UnionMT.__index = UnionMT

-- map(string,Struct) -> Union
function Union(choices)
  for name,fields in pairs(choices) do
    assert(type(name) == 'string')
    assert(isStruct(fields))
  end
  return setmetatable(UTIL.copyTable(choices), UnionMT)
end

-- SchemaT -> bool
local function isUnion(typ)
  return type(typ) == 'table' and getmetatable(typ) == UnionMT
end

-------------------------------------------------------------------------------

-- Struct = map(string,SchemaT)

-- SchemaT -> bool
function isStruct(typ)
  if type(typ) ~= 'table' then
    return false
  end
  for k,v in pairs(typ) do
    if type(k) ~= 'string' or not isSchemaT(v) then
      return false
    end
  end
  return true
end

-------------------------------------------------------------------------------

-- SchemaT = 'bool' | 'int' | 'double' | String
--         | Enum | Array | UpTo | Union | Struct

-- A -> bool
function isSchemaT(typ)
  return
    typ == bool or
    typ == int or
    typ == double or
    isString(typ) or
    isEnum(typ) or
    isArray(typ) or
    isUpTo(typ) or
    isUnion(typ) or
    isStruct(typ)
end

-------------------------------------------------------------------------------

-- SchemaT -> terralib.type
local function convertSchemaT(typ)
  if typ == bool then
    return bool
  elseif typ == int then
    return int
  elseif typ == double then
    return double
  elseif isString(typ) then
    return int8[typ.maxLen]
  elseif isEnum(typ) then
    return int
  elseif isArray(typ) then
    return convertSchemaT(typ.elemType)[typ.num]
  elseif isUpTo(typ) then
    return struct {
      length : uint32;
      values : convertSchemaT(typ.elemType)[typ.max];
    }
  elseif isUnion(typ) then
    local u = terralib.types.newstruct()
    u.entries:insert(terralib.newlist())
    for n,t in pairs(typ) do
      u.entries[1]:insert({field=n, type=convertSchemaT(t)})
    end
    local s = terralib.types.newstruct()
    s.entries:insert({field='type', type=int})
    s.entries:insert({field='u', type=u})
    return s
  elseif isStruct(typ) then
    local s = terralib.types.newstruct()
    for n,t in pairs(typ) do
      s.entries:insert({field=n, type=convertSchemaT(t)})
    end
    return s
  else assert(false) end
end

-- string, terralib.expr? -> terralib.quote
local function errorOut(msg, name)
  if name then
    return quote
      var stderr = C.fdopen(2, 'w')
      C.fprintf(stderr, '%s for option %s\n', msg, name)
      C.fflush(stderr)
      C.exit(1)
    end
  else
    return quote
      var stderr = C.fdopen(2, 'w')
      C.fprintf(stderr, '%s\n', msg)
      C.fflush(stderr)
      C.exit(1)
    end
  end
end

-- terralib.symbol, terralib.expr, terralib.expr, SchemaT -> terralib.quote
local function emitValueParser(name, lval, rval, typ)
  if typ == bool then
    return quote
      if [rval].type ~= JSON.json_boolean then
        [errorOut('Wrong type', name)]
      end
      [lval] = [bool]([rval].u.boolean)
    end
  elseif typ == int then
    return quote
      if [rval].type ~= JSON.json_integer then
        [errorOut('Wrong type', name)]
      end
      [lval] = [rval].u.integer
    end
  elseif typ == double then
    return quote
      if [rval].type ~= JSON.json_double then
        [errorOut('Wrong type', name)]
      end
      [lval] = [rval].u.dbl
    end
  elseif isString(typ) then
    return quote
      if [rval].type ~= JSON.json_string then
        [errorOut('Wrong type', name)]
      end
      if [rval].u.string.length >= [typ.maxLen] then
        [errorOut('String too long', name)]
      end
      C.strncpy([lval], [rval].u.string.ptr, [typ.maxLen])
    end
  elseif isEnum(typ) then
    return quote
      if [rval].type ~= JSON.json_string then
        [errorOut('Wrong type', name)]
      end
      var found = false
      escape for choice,value in pairs(typ) do emit quote
        if C.strcmp([rval].u.string.ptr, choice) == 0 then
          [lval] = value
          found = true
        end
      end end end
      if not found then
        [errorOut('Unexpected value', name)]
      end
    end
  elseif isArray(typ) then
    return quote
      if [rval].type ~= JSON.json_array then
        [errorOut('Wrong type', name)]
      end
      if [rval].u.array.length ~= [typ.num] then
        [errorOut('Wrong length', name)]
      end
      for i = 0,[typ.num] do
        var rval_i = [rval].u.array.values[i]
        [emitValueParser(name..'[i]', `[lval][i], rval_i, typ.elemType)]
      end
    end
  elseif isUpTo(typ) then
    return quote
      if [rval].type ~= JSON.json_array then
        [errorOut('Wrong type', name)]
      end
      if [rval].u.array.length > [typ.max] then
        [errorOut('Too many values', name)]
      end
      [lval].length = [rval].u.array.length
      for i = 0,[rval].u.array.length do
        var rval_i = [rval].u.array.values[i]
        [emitValueParser(name..'[i]', `[lval].values[i], rval_i, typ.elemType)]
      end
    end
  elseif isUnion(typ) then
    return quote
      if [rval].type ~= JSON.json_object then
        [errorOut('Wrong type', name)]
      end
      var foundType = false
      for i = 0,[rval].u.object.length do
        var nodeName = [rval].u.object.values[i].name
        var nodeValue = [rval].u.object.values[i].value
        if C.strcmp(nodeName, 'type') == 0 then
          foundType = true
          if nodeValue.type ~= JSON.json_string then
            [errorOut('Type field on union not a string', name)]
          end
          escape local j = 0; for choice,fields in pairs(typ) do emit quote
            if C.strcmp(nodeValue.u.string.ptr, [choice]) == 0 then
              [lval].type = [j]
              [emitValueParser(name, `[lval].u.[choice], rval, fields)]
              break
            end
          end j = j + 1; end end
          [errorOut('Unrecognized type on union', name)]
        end
      end
      if not foundType then
        [errorOut('Missing type field on union', name)]
      end
    end
  elseif isStruct(typ) then
    return quote
      var totalParsed = 0
      if [rval].type ~= JSON.json_object then
        [errorOut('Wrong type', name)]
      end
      for i = 0,[rval].u.object.length do
        var nodeName = [rval].u.object.values[i].name
        var nodeValue = [rval].u.object.values[i].value
        var parsed = false
        escape for fld,subTyp in pairs(typ) do emit quote
          if C.strcmp(nodeName, fld) == 0 then
            [emitValueParser(name..'.'..fld, `[lval].[fld], nodeValue, subTyp)]
            parsed = true
          end
        end end end
        if parsed then
          totalParsed = totalParsed + 1
        elseif C.strcmp(nodeName, 'type') ~= 0 then
          var stderr = C.fdopen(2, 'w')
          C.fprintf(
            stderr, ['Warning: Ignoring option '..name..'.%s\n'], nodeName)
        end
      end
      -- TODO: Assuming the json file contains no duplicate values
      if totalParsed < [UTIL.tableSize(typ)] then
        [errorOut('Missing options from config file')]
      end
    end
  else assert(false) end
end

-------------------------------------------------------------------------------

local ext = '.lua'
if #arg < 1 or not arg[1]:endswith(ext) then
  print('Usage: '..arg[0]..' <schema'..ext..'>')
  os.exit(1)
end
local baseName = arg[1]:sub(1, arg[1]:len() - ext:len())

local SCHEMA = dofile(arg[1])

local configStruct = convertSchemaT(SCHEMA.Config)
configStruct.name = 'Config'
local hdrFile = io.open(baseName..'.h', 'w')
hdrFile:write('// DO NOT EDIT THIS FILE, IT IS AUTOMATICALLY GENERATED\n')
hdrFile:write('\n')
hdrFile:write('#include <stdbool.h>\n')
hdrFile:write('#include <stdint.h>\n')
hdrFile:write('\n')
hdrFile:write(UTIL.prettyPrintStruct(configStruct, true)..';\n')
for name,typ in pairs(SCHEMA) do
  if isEnum(typ) then
    hdrFile:write('\n')
    hdrFile:write('typedef int '..name..';\n')
    for choice,value in pairs(typ) do
      hdrFile:write('#define '..name..'_'..choice..' '..value..'\n')
    end
  elseif isUnion(typ) then
    hdrFile:write('\n')
    hdrFile:write('typedef int '..name..';\n')
    local i = 0
    for choice,_ in pairs(typ) do
      hdrFile:write('#define '..name..'_'..choice..' '..i..'\n')
      i = i + 1
    end
  end
end
hdrFile:write('\n')
hdrFile:write('#ifdef __cplusplus\n')
hdrFile:write('extern "C" {\n')
hdrFile:write('#endif\n')
hdrFile:write('struct Config parse_config(const char* filename);\n')
hdrFile:write('#ifdef __cplusplus\n')
hdrFile:write('}\n')
hdrFile:write('#endif\n')
hdrFile:close()

local terra parse_config(filename : &int8) : configStruct
  var config : configStruct
  var f = C.fopen(filename, 'r');
  if f == nil then [errorOut('Cannot open config file')] end
  var res1 = C.fseek(f, 0, C.SEEK_END);
  if res1 ~= 0 then [errorOut('Cannot seek to end of config file')] end
  var len = C.ftell(f);
  if len < 0 then [errorOut('Cannot ftell config file')] end
  var res2 = C.fseek(f, 0, C.SEEK_SET);
  if res2 ~= 0 then [errorOut('Cannot seek to start of config file')] end
  var buf = [&int8](C.malloc(len));
  if buf == nil then [errorOut('Malloc error while parsing config')] end
  var res3 = C.fread(buf, 1, len, f);
  if res3 < len then [errorOut('Cannot read from config file')] end
  C.fclose(f);
  var errMsg : int8[256]
  var settings = JSON.json_settings{ 0, 0, nil, nil, nil, 0 }
  settings.settings = JSON.json_enable_comments
  var root = JSON.json_parse_ex(&settings, buf, len, errMsg)
  if root == nil then
    C.printf('JSON parsing error: %s\n', errMsg)
    C.exit(1)
  end
  [emitValueParser('Config', config, root, SCHEMA.Config)]
  JSON.json_value_free(root)
  C.free(buf)
  return config
end
terralib.saveobj(baseName..'.o', 'object', {parse_config=parse_config})
