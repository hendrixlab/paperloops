set xlabel "Angle, Rotation matrix to axis/angle"
set ylabel "Angle, Rotation matrix to quaternion to axis/angle"
set term post enhanced color solid "Helvetica 20"
set output "before_vs_after_angle.ps"
plot "test.txt" u 2:1 title "0-0, 0-1, and 1-0 junctions"

