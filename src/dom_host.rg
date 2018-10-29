-- Runs dom.rg standalone.
-- Reads configuration options in the same format as main simulation.
-- Uses default values for Ib and sigma.

-------------------------------------------------------------------------------
-- Imports
-------------------------------------------------------------------------------

import 'regent'

local DOM = require 'dom-desugared'
