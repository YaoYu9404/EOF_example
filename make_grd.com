#!/bin/csh

set n = 1
while ($n <= 12)
        gmt xyz2grd -R0/360/-81/81 -I30m vari.$n.xyz -Gvari.$n.grd
        set n = `awk 'BEGIN{print '$n'+1}'`
end

