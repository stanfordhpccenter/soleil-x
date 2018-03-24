local C = terralib.includecstring [[
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
]]
local JSON = terralib.includec('json.h')
local UTIL = require 'util'

-------------------------------------------------------------------------------

-- Enum = map(string,int)
-- SchemaT = Enum | bool | int | double | double[N] | map(string,SchemaT)

-- SchemaT -> bool
local function isEnum(typ)
  if type(typ) ~= 'table' then
    return false
  end
  for k,v in pairs(typ) do
    if type(k) ~= 'string' or type(v) ~= 'number' or v ~= math.floor(v) then
      return false
    end
  end
  return true
end

-- SchemaT -> bool
local function isDoubleVecT(typ)
  return terralib.types.istype(typ) and typ:isarray() and typ.type == double
end

-- A -> bool
local function isSchemaT(typ)
  if isEnum(typ) or
     typ == bool or
     typ == int or
     typ == double or
     isDoubleVecT(typ) then
    return true
  end
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

-- SchemaT -> terralib.type
local function convertSchemaT(typ)
  assert(isSchemaT(typ))
  if isEnum(typ) then
    return int
  elseif typ == bool then
    return bool
  elseif typ == int then
    return int
  elseif typ == double then
    return double
  elseif isDoubleVecT(typ) then
    return double[typ.N]
  else
    local s = terralib.types.newstruct()
    for n,t in pairs(typ) do
      s.entries:insert({field=n, type=convertSchemaT(t)})
    end
    return s
  end
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
  assert(isSchemaT(typ))
  if isEnum(typ) then
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
  elseif typ == bool then
    return quote
      if [rval].type ~= JSON.json_boolean then
        [errorOut('Wrong type', name)]
      end
      [lval] = [rval].u.boolean
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
  elseif isDoubleVecT(typ) then
    return quote
      if [rval].type ~= JSON.json_array then
        [errorOut('Wrong type', name)]
      end
      if [rval].u.array.length ~= [typ.N] then
        [errorOut('Wrong length', name)]
      end
      for i = 0,[typ.N] do
        var rval_i = [rval].u.array.values[i]
        if rval_i.type ~= JSON.json_double then
          [errorOut('Wrong element type', name)]
        end
        [lval][i] = rval_i.u.dbl
      end
    end
  else
    return quote
      var totalParsed = 0
      if [rval].type ~= JSON.json_object then
        [errorOut('Wrong type', name)]
      end
      for i = 0,[rval].u.object.length do
        var nodeName = [rval].u.object.values[i].name
        var nodeValue = [rval].u.object.values[i].value
        var parsed = false
        escape for fld,subType in pairs(typ) do emit quote
          if C.strcmp(nodeName, fld) == 0 then
            [emitValueParser(nodeName, `[lval].[fld], nodeValue, subType)]
            parsed = true
          end
        end end end
        if parsed then
          totalParsed = totalParsed + 1
        else
          var stderr = C.fdopen(2, 'w')
          C.fprintf(stderr, 'Warning: Ignoring option %s\n', nodeName)
        end
      end
      -- TODO: Assuming the json file contains no duplicate values
      if totalParsed < [UTIL.tableSize(typ)] then
        [errorOut('Missing options from config file')]
      end
    end
  end
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
hdrFile:write(UTIL.prettyPrintStruct(configStruct, true)..';\n')
for name,typ in pairs(SCHEMA) do
  if isEnum(typ) then
    hdrFile:write('\n')
    hdrFile:write('typedef int '..name..';\n')
    for choice,value in pairs(typ) do
      hdrFile:write('#define '..name..'_'..choice..' '..value..'\n')
    end
  end
end
hdrFile:write('\n')
hdrFile:write('#ifdef __cplusplus\n')
hdrFile:write('extern "C" {\n')
hdrFile:write('#endif\n')
hdrFile:write('struct Config parse_config(char* filename);\n')
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
  [emitValueParser('<root>', config, root, SCHEMA.Config)]
  JSON.json_value_free(root)
  C.free(buf)
  return config
end
terralib.saveobj(baseName..'.o', 'object', {parse_config=parse_config})
