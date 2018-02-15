-- Load a config schema header as a module, and augment it with '==' operators
-- on enum-type fields, as well as a 'parse' function that accepts compatible
-- config files in JSON format.
-- NOTE:
-- * The schema header must only declare a single struct, named 'Config'
--   (other structs must be nested inside it).
-- * The schema header must use the MK_ENUM macro from 'enum.h' to declare
--   enums.

-------------------------------------------------------------------------------

import 'regent'

local C = terralib.includecstring [[
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
]]
local JSON = terralib.includec('json.h')
local UTIL = require 'util'

local Exports = {}

-------------------------------------------------------------------------------

-- A -> bool
local function isEnum(x)
  if not terralib.types.istype(x) or not x:isstruct() then
    return false
  end
  for _,e in ipairs(x.entries) do
    local n,t = UTIL.parseStructEntry(e)
    if n ~= '__value' or t ~= int then
      return false
    end
  end
  return true
end

-- string, terralib.expr? -> terralib.quote
local function errorOut(msg, name)
  if name then
    return quote
      C.printf('%s for option %s\n', msg, name)
      C.exit(1)
    end
  else
    return quote
      C.printf('%s\n', msg)
      C.exit(1)
    end
  end
end

-- terralib.symbol, terralib.expr, terralib.expr, terralib.type
--   -> terralib.quote
local function emitValueParser(name, lval, rval, type)
  if isEnum(type) then
    return quote
      if [rval].type ~= JSON.json_string then
        [errorOut('Wrong type', name)]
      end
      var found = false
      escape for i,choice in ipairs(type.__choices) do emit quote
        if C.strcmp([rval].u.string.ptr, choice) == 0 then
          [lval].__value = i-1
          found = true
        end
      end end end
      if not found then
        [errorOut('Unexpected value', name)]
      end
    end
  elseif type == int then
    return quote
      if [rval].type ~= JSON.json_integer then
        [errorOut('Wrong type', name)]
      end
      [lval] = [rval].u.integer
    end
  elseif type == double then
    return quote
      if [rval].type ~= JSON.json_double then
        [errorOut('Wrong type', name)]
      end
      [lval] = [rval].u.dbl
    end
  elseif type:isarray() and type.type == double then
    return quote
      if [rval].type ~= JSON.json_array then
        [errorOut('Wrong type', name)]
      end
      if [rval].u.array.length ~= [type.N] then
        [errorOut('Wrong length', name)]
      end
      for i = 0,[type.N] do
        var rval_i = [rval].u.array.values[i]
        if rval_i.type ~= JSON.json_double then
          [errorOut('Wrong element type', name)]
        end
        [lval][i] = rval_i.u.dbl
      end
    end
  elseif type:isstruct() then
    local entries = {}
    local numEntries = 0
    for _,e in ipairs(type.entries) do
      local n,t = UTIL.parseStructEntry(e)
      entries[n] = t
      numEntries = numEntries + 1
    end
    return quote
      var totalParsed = 0
      if [rval].type ~= JSON.json_object then
        [errorOut('Wrong type', name)]
      end
      for i = 0,[rval].u.object.length do
        var nodeName = [rval].u.object.values[i].name
        var nodeValue = [rval].u.object.values[i].value
        var parsed = false
        escape for fld,subType in pairs(entries) do emit quote
          if C.strcmp(nodeName, fld) == 0 then
            [emitValueParser(nodeName, `[lval].[fld], nodeValue, subType)]
            parsed = true
          end
        end end end
        if parsed then
          totalParsed = totalParsed + 1
        else
          C.printf('Ignoring option %s\n', nodeName)
        end
      end
      -- TODO: Assuming the json file contains no duplicate values
      if totalParsed < [numEntries] then
        [errorOut('Missing options from config file')]
      end
    end
  else assert(false) end
end

-- terralib.struct -> terralib.function
local function emitConfigParser(configStruct)
  local terra parseConfig(filename : &int8) : configStruct
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
    [emitValueParser('<root>', config, root, configStruct)]
    JSON.json_value_free(root)
    C.free(buf)
    return config
  end
  return parseConfig
end

-- string -> terralib.module
function Exports.processSchema(schemaHdr)
  local configMod = terralib.includec(schemaHdr)
  for enumName,enumT in pairs(configMod) do
    if isEnum(enumT) then
      enumT.__choices = {}
      for key,val in pairs(configMod) do
        if key:startswith(enumName..'_') then
          local choice = key:sub(enumName:len()+2)
          enumT.__choices[val+1] = choice
          enumT[choice] = val
        end
      end
      enumT.metamethods.__eq = macro(function(a,b) return `a.__value == b end)
    end
  end
  configMod.parseConfig = emitConfigParser(configMod.Config)
  return configMod
end

-------------------------------------------------------------------------------

return Exports
