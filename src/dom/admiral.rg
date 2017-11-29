-- Proxies the functionality of admiral.t, for running dom.rg standalone.

import 'regent'

local A = {}

A.primColors = terralib.memoize(function()
  return regentlib.newsymbol()
end)

A.configSymbol = terralib.memoize(function()
  return regentlib.newsymbol()
end)

function A.registerTask(tsk, name) end

function A.registerFun(fun, name) end

function A.registerStruct(s) end

return A
