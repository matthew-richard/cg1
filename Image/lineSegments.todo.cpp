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
	float numerator = (y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1;
	return numerator / length();
}
void  OrientedLineSegment::getPerpendicular(float& x,float &y) const
{
	float px = (y2 - y1);
	float py = (x1 - x2);
	float norm = sqrt(px * px + py * py);
	x = px / norm;
	y = py / norm;
}

void  OrientedLineSegment::GetSourcePosition(const OrientedLineSegment& source,const OrientedLineSegment& destination,
											 const int& targetX,const int& targetY,
											 float& sourceX,float& sourceY)
{
	sourceX=sourceY=0;
}
