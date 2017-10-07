-- Proxies the functionality of admiral.t, for running dom.rg standalone.

-------------------------------------------------------------------------------
-- Imports
-------------------------------------------------------------------------------

import 'regent'

local A = {}

-------------------------------------------------------------------------------
-- Proxy Admiral functionality
-------------------------------------------------------------------------------

A.primColors = terralib.memoize(function()
  return regentlib.newsymbol()
end)

function A.primPartDims()
  local dop = regentlib.config['parallelize-dop']
  local NX,NY,NZ = dop:match('^(%d+),(%d+),(%d+)$')
  if NX then
    return tonumber(NX),tonumber(NY),tonumber(NZ)
  end
  local num = assert(tonumber(dop))
  local factors = terralib.newlist()
  while num > 1 do
    for p = 2, num do
      if num % p == 0 then
        factors:insert(p)
        num = num / p
        break
      end
    end
  end
  NX,NY,NZ = 1,1,1
  for i = 1, #factors do
    if i % 3 == 1 then NX = NX * factors[i] end
    if i % 3 == 2 then NY = NY * factors[i] end
    if i % 3 == 0 then NZ = NZ * factors[i] end
  end
  return NX,NY,NZ
end

return A
