function d = point_to_line(pt, v1, v2)
% Computes the distance from a point to a line
      a = v1 - v2;
      b = pt - v2;
      d = norm(cross(a,b)) / norm(a);