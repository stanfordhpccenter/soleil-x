-- Runs dom.rg standalone.
-- Reads configuration options in the same format as main simulation.
-- Uses default values for Ib and sigma.

-------------------------------------------------------------------------------
-- Imports
-------------------------------------------------------------------------------

import 'regent'

local SCHEMA = terralib.includec("config_schema.h")

local MAX_ANGLES_PER_QUAD = 44

-------------------------------------------------------------------------------
-- Proxy radiation grid
-------------------------------------------------------------------------------

struct Point_columns {
  G : double;
  S : double;
  Ib : double;
  sigma : double;
}

-------------------------------------------------------------------------------
-- Import DOM module
-------------------------------------------------------------------------------

local DOM = (require 'dom-desugared')(MAX_ANGLES_PER_QUAD, Point_columns, SCHEMA)
