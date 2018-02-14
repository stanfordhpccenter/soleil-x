-- Generate code for dumping/loading a subset of fields to/from an HDF file.
-- NOTE:
-- * Both functions require an intermediate region to perform the data
--   transfer. This region 's' must have the same size as 'r', and must be
--   partitioned in the same way.
-- * The dimensions will be flipped in the output file.
-- * You need to link to the HDF library to use these functions.

-------------------------------------------------------------------------------

import 'regent'

local Exports = {}

local UTIL = require 'util'

local HDF5 = terralib.includec(assert(os.getenv('HDF_HEADER')))
-- HACK: Hardcoding missing #define's
HDF5.H5F_ACC_TRUNC = 2
HDF5.H5P_DEFAULT = 0

-------------------------------------------------------------------------------
-- Dumping
-------------------------------------------------------------------------------

-- regentlib.index_type, regentlib.index_type, terralib.struct, string*
--   -> regentlib.task
function Exports.mkDump(indexType, colorType, fSpace, flds)
  flds = terralib.newlist(flds)

  local terra create(fname : rawstring, hiBound : indexType)
    var fid = HDF5.H5Fcreate(fname, HDF5.H5F_ACC_TRUNC,
                             HDF5.H5P_DEFAULT, HDF5.H5P_DEFAULT)
    var dataSpace : int32
    escape
      if indexType == int1d then
        emit quote
          var sizes : HDF5.hsize_t[1]
          sizes[0] = hiBound.__ptr + 1
          dataSpace = HDF5.H5Screate_simple(1, sizes, [&uint64](0))
        end
      elseif indexType == int2d then
        emit quote
          -- Legion defaults to column-major layout, so we have to reverse.
          var sizes : HDF5.hsize_t[2]
          sizes[1] = hiBound.__ptr.x + 1
          sizes[0] = hiBound.__ptr.y + 1
          dataSpace = HDF5.H5Screate_simple(2, sizes, [&uint64](0))
        end
      elseif indexType == int3d then
        emit quote
          -- Legion defaults to column-major layout, so we have to reverse.
          var sizes : HDF5.hsize_t[3]
          sizes[2] = hiBound.__ptr.x + 1
          sizes[1] = hiBound.__ptr.y + 1
          sizes[0] = hiBound.__ptr.z + 1
          dataSpace = HDF5.H5Screate_simple(3, sizes, [&uint64](0))
        end
      else assert(false) end
      local header = terralib.newlist() -- terralib.quote*
      local footer = terralib.newlist() -- terralib.quote*
      -- terralib.type -> terralib.quote
      local function toHType(T)
        -- TODO: Not supporting: pointers, vectors, non-primitive arrays
        if T:isprimitive() then
          return
            -- HACK: Hardcoding missing #define's
            (T == int)    and HDF5.H5T_STD_I32LE_g  or
            (T == int8)   and HDF5.H5T_STD_I8LE_g   or
            (T == int16)  and HDF5.H5T_STD_I16LE_g  or
            (T == int32)  and HDF5.H5T_STD_I32LE_g  or
            (T == int64)  and HDF5.H5T_STD_I64LE_g  or
            (T == uint)   and HDF5.H5T_STD_U32LE_g  or
            (T == uint8)  and HDF5.H5T_STD_U8LE_g   or
            (T == uint16) and HDF5.H5T_STD_U16LE_g  or
            (T == uint32) and HDF5.H5T_STD_U32LE_g  or
            (T == uint64) and HDF5.H5T_STD_U64LE_g  or
            (T == bool)   and HDF5.H5T_STD_U8LE_g   or
            (T == float)  and HDF5.H5T_IEEE_F32LE_g or
            (T == double) and HDF5.H5T_IEEE_F64LE_g or
            assert(false)
        elseif T:isarray() then
          local elemType = toHType(T.type)
          local arrayType = symbol(HDF5.hid_t, 'arrayType')
          header:insert(quote
            var dims : HDF5.hsize_t[1]
            dims[0] = T.N
            var elemType = [elemType]
            var [arrayType] = HDF5.H5Tarray_create2(elemType, 1, dims)
          end)
          footer:insert(quote
            HDF5.H5Tclose(arrayType)
          end)
          return arrayType
        else assert(false) end
      end
      -- terralib.struct, set(string), string -> ()
      local function emitFieldDecls(fs, whitelist, prefix)
        -- TODO: Only supporting pure structs, not fspaces
        assert(fs:isstruct())
        for _,e in ipairs(fs.entries) do
          local name, type = UTIL.parseStructEntry(e)
          if whitelist and not whitelist[name] then
            -- do nothing
          elseif type == int2d then
            -- Hardcode special case: int2d structs are stored packed
            local hName = prefix..name
            local int2dType = symbol(HDF5.hid_t, 'int2dType')
            local dataSet = symbol(HDF5.hid_t, 'dataSet')
            header:insert(quote
              var [int2dType] = HDF5.H5Tcreate(HDF5.H5T_COMPOUND, 16)
              HDF5.H5Tinsert(int2dType, "x", 0, HDF5.H5T_STD_I64LE_g)
              HDF5.H5Tinsert(int2dType, "y", 8, HDF5.H5T_STD_I64LE_g)
              var [dataSet] = HDF5.H5Dcreate2(
                fid, hName, int2dType, dataSpace,
                HDF5.H5P_DEFAULT, HDF5.H5P_DEFAULT, HDF5.H5P_DEFAULT)
            end)
            footer:insert(quote
              HDF5.H5Dclose(dataSet)
              HDF5.H5Tclose(int2dType)
            end)
          elseif type == int3d then
            -- Hardcode special case: int3d structs are stored packed
            local hName = prefix..name
            local int3dType = symbol(HDF5.hid_t, 'int3dType')
            local dataSet = symbol(HDF5.hid_t, 'dataSet')
            header:insert(quote
              var [int3dType] = HDF5.H5Tcreate(HDF5.H5T_COMPOUND, 24)
              HDF5.H5Tinsert(int3dType, "x", 0, HDF5.H5T_STD_I64LE_g)
              HDF5.H5Tinsert(int3dType, "y", 8, HDF5.H5T_STD_I64LE_g)
              HDF5.H5Tinsert(int3dType, "z", 16, HDF5.H5T_STD_I64LE_g)
              var [dataSet] = HDF5.H5Dcreate2(
                fid, hName, int3dType, dataSpace,
                HDF5.H5P_DEFAULT, HDF5.H5P_DEFAULT, HDF5.H5P_DEFAULT)
            end)
            footer:insert(quote
              HDF5.H5Dclose(dataSet)
              HDF5.H5Tclose(int3dType)
            end)
          elseif type:isstruct() then
            emitFieldDecls(type, nil, prefix..name..'.')
          else
            local hName = prefix..name
            local hType = toHType(type)
            local dataSet = symbol(HDF5.hid_t, 'dataSet')
            header:insert(quote
              var hType = [hType]
              var [dataSet] = HDF5.H5Dcreate2(
                fid, hName, hType, dataSpace,
                HDF5.H5P_DEFAULT, HDF5.H5P_DEFAULT, HDF5.H5P_DEFAULT)
            end)
            footer:insert(quote
              HDF5.H5Dclose(dataSet)
            end)
          end
        end
      end
      emitFieldDecls(fSpace, flds:toSet(), '')
      emit quote [header] end
      emit quote [footer:reverse()] end
    end
    HDF5.H5Sclose(dataSpace)
    HDF5.H5Fclose(fid)
  end

  local task dump(colors : ispace(colorType),
                  filename : int8[256],
                  r : region(ispace(indexType),fSpace),
                  s : region(ispace(indexType),fSpace),
                  p_r : partition(disjoint, r, colors),
                  p_s : partition(disjoint, s, colors))
  where reads(r.[flds]), reads writes(s.[flds]), r * s do
    -- TODO: Sanity checks: bounds.lo == 0, same size
    create(filename, r.bounds.hi)
    attach(hdf5, s.[flds], filename, regentlib.file_read_write)
    for c in colors do
      var p_r_c = p_r[c]
      var p_s_c = p_s[c]
      acquire(p_s_c.[flds])
      copy(p_r_c.[flds], p_s_c.[flds])
      release(p_s_c.[flds])
    end
    detach(hdf5, s.[flds])
  end

  return dump
end

-------------------------------------------------------------------------------
-- Loading
-------------------------------------------------------------------------------

-- regentlib.index_type, regentlib.index_type, terralib.struct, string*
--   -> regentlib.task
function Exports.mkLoad(indexType, colorType, fSpace, flds)
  flds = terralib.newlist(flds)

  local task load(colors : ispace(colorType),
                  filename : int8[256],
                  r : region(ispace(indexType),fSpace),
                  s : region(ispace(indexType),fSpace),
                  p_r : partition(disjoint, r, colors),
                  p_s : partition(disjoint, s, colors))
  where reads writes(r.[flds]), reads writes(s.[flds]), r * s do
    -- TODO: Sanity checks: bounds.lo == 0, same size
    attach(hdf5, s.[flds], filename, regentlib.file_read_only)
    for c in colors do
      var p_r_c = p_r[c]
      var p_s_c = p_s[c]
      acquire(p_s_c.[flds])
      copy(p_s_c.[flds], p_r_c.[flds])
      release(p_s_c.[flds])
    end
    detach(hdf5, s.[flds])
  end

  return load
end

-------------------------------------------------------------------------------

return Exports
