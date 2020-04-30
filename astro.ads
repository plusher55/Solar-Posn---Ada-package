package Astro is

--  ----------------------------------------------------------------------------
--  Author: Paul W Lusher   winter 2016/2017
--
--    sunPosition calculates the local azimuth and elevation of the
--    sun at a specific location and time using equations from the
--    US Naval Observatory sites:
--      http://aa.usno.navy.mil/faq/docs/SunApprox.php
--      http://aa.usno.navy.mil/faq/docs/GAST.php
--      http://aa.usno.navy.mil/faq/docs/Alt_Az.php
--
--    This provides an approximation to equations from:
--    The Astronomical Almanac for the Year 2010.
--    U.S. Govt. Printing Office. p. C5. ISBN 978-0-7077-4082-9.
--
--  ----------------------------------------------------------------------------
--      input parameters
--        year =year, e.g., 2015
--        day  =day of year (1-366, e.g. feb 1=32)
--        hour =hours plus fraction in UTC, e.g. 1430=14.50
--        lat  =latitude in degrees (north is positive)
--        long =longitude in degrees (east is positive)
--
--      output parameters all degrees except solar distance
--        Az   = sun azimuth angle (measured east from north, 0 to 360 degs)
--        Elev = sun elevation angle (degs)
--        hourAngle   = solar hour angle
--        Declination = declination
--        rAscension  = right ascension
--
--  ----------------------------------------------------------------------------
--
--  Results have been verified by checking various locations with differnt
--  time values (specific or day-long).
--  Verification done by various web-based solar calculators, including:
--
--       http://aa.usno.navy.mil/data/docs/AltAz.php
--       https://www.suncalc.org
--       http://www.satellite-calculations.com/Satellite/suncalc.htm
--
--  ----------------------------------------------------------------------------
procedure sunPosition(
 year : in integer;
 day  : in integer;
 hour : in float;
 lat  : in float;
 long : in float;
 ----------
 Az            : out float;
 Elev          : out float;
 hourAngle     : out float;
 Declination   : out float;
 rAscension    : out float );

end Astro;
