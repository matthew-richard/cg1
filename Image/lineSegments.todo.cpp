#include "lineSegments.h"
#include <math.h>
#include <iostream>

////////////////////////////
// Image processing stuff //
////////////////////////////
float OrientedLineSegment::length(void) const
{
	return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
}
float OrientedLineSegment::distance(const int& x,const int& y) const
{
	float px, py;
	getPerpendicular(px, py);
	float line_dist = dot(x - x1, y - y1, px, py);

	float dot1 = dot(x - x1, y - y1, x2 - x1, y2 - y1);
	float dot2 = dot(x - x2, y - y2, x1 - x2, y1 - y2);
	if (dot1 < 0)
	{
		return dist(x1, y1, x, y);
	}
	else if (dot2 < 0)
	{
		return dist(x2, y2, x, y);
	}
	else
	{
		return line_dist;
	}
}
void  OrientedLineSegment::getPerpendicular(float& x,float &y) const
{
	float px = (y2 - y1);
	float py = (x1 - x2);
	float n = norm(px, py);
	x = px / n;
	y = py / n;
}

void  OrientedLineSegment::GetSourcePosition(const OrientedLineSegment& source,const OrientedLineSegment& destination,
											 const int& targetX,const int& targetY,
											 float& sourceX,float& sourceY)
{
	float px, py; destination.getPerpendicular(px, py);

	// Signed perpendicular distance from target point to destination line
	float line_dist = dot(targetX - destination.x1, targetY - destination.y1, px, py);

	// Target point's position along line segment. 0 means the point is at destination's first point, 1 means it's at destination's second
	// point, outside of [0, 1] means point lies "outside" of line segment
	float line_pos = dot(destination.x2 - destination.x1, destination.y2 - destination.y1, targetX - destination.x1, targetY - destination.y1)
					   / pow(destination.length(), 2);

	float px_s, py_s; source.getPerpendicular(px_s, py_s);
	sourceX = source.x1 + line_pos * (source.x2 - source.x1) + line_dist * px_s;
	sourceY = source.y1 + line_pos * (source.y2 - source.y1) + line_dist * py_s;
}

float dot(const float& x1, const float& y1, const float& x2, const float& y2)
{
	return x1 * x2 + y1 * y2;
}

float dist(const float& x1, const float& y1, const float& x2, const float& y2)
{
	return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

float norm(const float& x, const float& y)
{
	sqrt(x * x + y * y);
}
