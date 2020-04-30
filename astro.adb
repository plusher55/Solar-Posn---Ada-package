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
--  ----------------------------------------------------------------------------
-- From: https://aa.usno.navy.mil/faq/docs/SunApprox.php:
--
-- "Given below is a simple algorithm for computing the Sun's angular
-- coordinates to an accuracy of about 1 arcminute within two centuries of 2000.
-- The algorithm's accuracy degrades gradually beyond its four-century window
-- of applicability. This accuracy is quite adequate for computing, for example,
-- the times of sunrise and sunset, or solar transit. For navigational
-- purposes it would provide about 1 nautical mile accuracy."
--
--  ----------------------------------------------------------------------------

with Ada.Text_IO;             use Ada.Text_IO;
with Ada.Float_Text_IO;       use Ada.Float_Text_IO;
with Ada.Long_Float_Text_IO;  use Ada.Long_Float_Text_IO;
with Ada.Numerics.Elementary_Functions;
     use Ada.Numerics.Elementary_Functions;
with Ada.Exceptions;              use  Ada.Exceptions;

package body Astro is

-------------------------------------------------------------
function Limit_Degrees_360 (Theta : float) return float is
   Limited_Theta : float := Theta;
begin
   if Limited_Theta >= 360.0 then
      loop
         Limited_Theta := Limited_Theta - 360.0;
         if (Limited_Theta < 360.0) then
            exit;
         end if;
      end loop;

   elsif Limited_Theta < 0.0 then
      loop
         Limited_Theta := Limited_Theta + 360.0;
         if (Limited_Theta >= 0.0) then
            exit;
         end if;
      end loop;
   end if;

   return Limited_Theta;
end Limit_Degrees_360;


--  ----------------------------------------------------------------------------
-- UTC 2013, 1 Jan, 00:30 generates 2,456,293.520833
--  matches https://en.wikipedia.org/wiki/Julian_day example
--
function Days_From_J2000p0 (year : integer;
                            day  : integer;
                            hour : float) return Float is

   delta2000, leap : integer := 0;
   JD2000          : Long_Float := 2_451_545.0;

begin
   delta2000 := year-2000;
   -- valid 2001 thru 2099 b/c of leap yr logic
   if delta2000 > 0 then
      leap    := (delta2000-1)/4 + 1;
   end if;

   return float(delta2000)*365.0 - 0.5 + float(leap) + float(day-1) + (hour/24.0);

end Days_From_J2000p0;


--  ----------------------------------------------------------------------------
function Previous_Midnight (days  : Float) return Float is

   days_int, days_frc : float;

begin
   -- get previous midnight
   days_int := Float(Integer(days));
   days_frc := days - days_int;
   if days_frc < 0.5 then
      days_int := days_int - 1.0;
   end if;
   return  days_int + 0.5;

end Previous_Midnight;



--  ----------------------------------------------------------------------------
--  Author: Paul W Lusher   winter 2016/2017
--
--    This subroutine calculates the local azimuth and elevation of the sun at
--    a specific location and time.
--
--    The UTC to Julian conversion is only correct for years 2001 - 2099
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
 rAscension    : out float ) is

   pi        : constant float := Ada.Numerics.PI;
   twopi     : constant float := 2.0 * pi;
   deg_rad   : constant float := pi/180.0;

   days, days0     : float;

   mLongitude, mAnomaly, eLongitude, OblqEcliptic : float;
   raNumerator, raDenominator : float;
   gmtSidereal, lmtSidereal : float;
   latrad, refrac : float;

begin
-- change latitude to radians
   latrad := lat*deg_rad;

-- time in days since NOON Jan 1 2000, equation valid 2001 thru 2099
   days   := Days_From_J2000p0 (year, day, hour);
   days0  := Previous_Midnight (days);  -- back to previous midnight
----------------------------

-- mean anomaly
   mAnomaly := 357.529 + 0.985_600_28*days;
   mAnomaly := Limit_Degrees_360(mAnomaly);

-- mean longitude
   mLongitude := 280.459 + 0.985_6474*days;
   mLongitude := Limit_Degrees_360(mLongitude);

-- calculate ecliptic longitude & obliquity of ecliptic (radians)
   mAnomaly   := mAnomaly*deg_rad;  -- need radians between 0 and 2*pi
   eLongitude := mLongitude + 1.915*sin(mAnomaly) + 0.020*sin(2.0*mAnomaly);
   eLongitude := Limit_Degrees_360(eLongitude);
   eLongitude := eLongitude*deg_rad;

   OblqEcliptic := 23.439-0.000_000_36*days; -- mean obliquity of the ecliptic
   OblqEcliptic := OblqEcliptic*deg_rad;     -- radians 0 to 2*pi

-- calculate right ascension & declination (radians)
   raNumerator := cos(OblqEcliptic)*sin(eLongitude);
   raDenominator := cos(eLongitude);
   rAscension  := arctan(raNumerator, raDenominator);

   Declination := arcsin(sin(OblqEcliptic)*sin(eLongitude));

-- calculate Greenwich mean sidereal time in hours
   gmtSidereal := 6.697_374_558 + 0.065_709_824_419*days0
                  + 1.002_737_909_35 * hour;
   gmtSidereal := Float'Remainder(gmtSidereal,24.0);
   if gmtSidereal < 0.0 then gmtSidereal := gmtSidereal+24.0; end if;
-- *** Note: use mean vs. apparent sidereal time, introducing
--           an error of about 1 second.

-- calculate local mean sidereal time (radians)
   lmtSidereal := gmtSidereal + long/15.0;
   lmtSidereal := Float'Remainder(lmtSidereal,24.0);
   if lmtSidereal < 0.0 then lmtSidereal := lmtSidereal+24.0; end if;

   lmtSidereal := lmtSidereal*15.0*deg_rad;

-- calculate hour angle in radians -pi to pi
   hourAngle := Float(lmtSidereal)-rAscension;
   if hourAngle < -pi then hourAngle := hourAngle+twopi; end if;
   if hourAngle >  pi then hourAngle := hourAngle-twopi; end if;

-- calculate azimuth & elevation
   Elev := arcsin(sin(Declination)*sin(latrad)+cos(Declination)*cos(latrad)*cos(hourAngle));
   Az := arcsin(-cos(Declination)*sin(hourAngle)/cos(Elev));

   if sin(Declination)-sin(Elev)*sin(latrad) < 0.0 then
      Az := pi-Az;
   end if;

   Elev := Elev/deg_rad;   -- to degrees for final output

   -- account for refraction (remove this code if apparent elev not needed)
   --   below, constant 3.51823=1013.25 mb/288 C (US Standard atmosphere)
   if(Elev >= 19.225) then
      refrac := 0.00452*3.51823/tan(Elev*deg_rad);
   elsif (Elev > -0.766 and Elev < 19.225) then
      refrac := 3.51823*(0.1594+0.0196*Elev+0.00002*Elev**2)/(1.0+0.505*Elev+0.0845*Elev**2);
   elsif (Elev <= -0.766) then
      refrac := 0.0;
   end if;
   Elev := Elev+refrac;

-- return degrees
   Az := Limit_Degrees_360(Az/deg_rad);
   hourAngle   := hourAngle/deg_rad;
   Declination := Declination/deg_rad;
   rAscension  := Limit_Degrees_360(rAscension/deg_rad);

   exception
      when exp: others =>
         Put("Solar Posn exception " & Exception_Name( exp ) );
         New_Line;
end sunPosition;

end Astro;
