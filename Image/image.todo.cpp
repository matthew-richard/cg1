#include "image.h"
#include <stdlib.h>
#include <math.h>

////////////////////////////
// Image processing stuff //
////////////////////////////
double Gaussian(const float& variance, const double& x);
float Distance(const float& x1, const float& y1, const float& x2, const float& y2);
float Clamp32(const float& x);
float ClampF(const float& x);
float Mirror(const int& x, const int& cutoff);
void SetRotatedSize(const Image32 * inputImage, const float& angle, Image32& outputImage, float& centerX, float& centerY);


Pixel::Pixel(const Pixel32& p)
{
	a = p.a / (float)255;
	r = p.r / (float)255;
	g = p.g / (float)255;
	b = p.b / (float)255;
}
Pixel32::Pixel32(const Pixel& p)
{
	a = (int)round(p.a * 255);
	r = (int)round(p.r * 255);
	g = (int)round(p.g * 255);
	b = (int)round(p.b * 255);
}

int Image32::AddRandomNoise(const float& noise,Image32& outputImage) const
{
	outputImage.setSize(width(), height());
	for (int u = 0; u < width(); u++)
	for (int v = 0; v < height(); v++)
	{
		outputImage.pixel(u, v).a = 255; // Not applying noise to alpha

		// Add random number between [-noise * 255, noise * 255]
		int rand_r = pixel(u, v).r + ((rand() / (float)RAND_MAX) - 0.5) * 2 * (noise * 255);
		int rand_g = pixel(u, v).g + ((rand() / (float)RAND_MAX) - 0.5) * 2 * (noise * 255);
		int rand_b = pixel(u, v).b + ((rand() / (float)RAND_MAX) - 0.5) * 2 * (noise * 255);
		outputImage.pixel(u, v).r = Clamp32(rand_r);
		outputImage.pixel(u, v).g = Clamp32(rand_g);
		outputImage.pixel(u, v).b = Clamp32(rand_b);
	}
	return 1;
}
int Image32::Brighten(const float& brightness,Image32& outputImage) const
{
	outputImage.setSize(w, h);
	for (int y = 0; y < h; y++) {
		for (int x = 0; x < w; x++)
		{
			const Pixel32& pixelIn = this->pixel(x, y);
			Pixel32& pixelOut = outputImage(x, y);
			pixelOut.a = pixelIn.a;

			// scale rgb and clamp to [0, 255]
			pixelOut.r = Clamp32(pixelIn.r * brightness);
			pixelOut.g = Clamp32(pixelIn.g * brightness);
			pixelOut.b = Clamp32(pixelIn.b * brightness);
		}
	}
	return 1;
}

int Image32::Luminance(Image32& outputImage) const
{
	outputImage.setSize(w, h);
	for (int y = 0; y < h; y++) {
		for (int x = 0; x < w; x++)
		{
			const Pixel32& pixelIn = this->pixel(x, y);
			unsigned char luminance = round(0.3 * pixelIn.r + 0.59 * pixelIn.g + 0.11 * pixelIn.b);

			Pixel32& pixelOut = outputImage(x, y);
			pixelOut.a = pixelIn.a;
			pixelOut.r = pixelOut.g = pixelOut.b = luminance;
		}
	}
	return 1;
}

int Image32::Contrast(const float& contrast,Image32& outputImage) const
{
	outputImage.setSize(w, h);

	// Compute average luminance
	float avgLuminance = 0;
	for (int y = 0; y < h; y++) {
		for (int x = 0; x < w; x++)
		{
			const Pixel32& pixelIn = this->pixel(x, y);
			float luminance = 0.3 * pixelIn.r + 0.59 * pixelIn.g + 0.11 * pixelIn.b;
			avgLuminance += luminance / (w * h);
		}
	}

	for (int y = 0; y < h; y++) {
		for (int x = 0; x < w; x++)
		{
			const Pixel32& pixelIn = this->pixel(x, y);
			Pixel32& pixelOut = outputImage(x, y);
			pixelOut.a = pixelIn.a;

			// Adjust contrast and clamp to [0, 255]
			pixelOut.r = Clamp32(avgLuminance + contrast * (pixelIn.r - avgLuminance));
			pixelOut.g = Clamp32(avgLuminance + contrast * (pixelIn.g - avgLuminance));
			pixelOut.b = Clamp32(avgLuminance + contrast * (pixelIn.b - avgLuminance));
		}
	}
	return 1;
}

int Image32::Saturate(const float& saturation,Image32& outputImage) const
{
	outputImage.setSize(width(), height());
	for (int u = 0; u < width(); u++)
	for (int v = 0; v < height(); v++)
	{
		const Pixel32& p = pixel(u, v);
		Pixel32& o = outputImage.pixel(u, v);

		float luminance = 0.3 * p.r + 0.59 * p.g + 0.11 * p.b;

		o.a = p.a;
		o.r = Clamp32(luminance + (p.r - luminance) * saturation);
		o.g = Clamp32(luminance + (p.g - luminance) * saturation);
		o.b = Clamp32(luminance + (p.b - luminance) * saturation);
	}
	
	return 1;
}

int Image32::Quantize(const int& bits,Image32& outputImage) const
{
	return 0;
}

int Image32::RandomDither(const int& bits,Image32& outputImage) const
{
	return 0;
}
int Image32::OrderedDither2X2(const int& bits,Image32& outputImage) const
{
	return 0;
}

int Image32::FloydSteinbergDither(const int& bits,Image32& outputImage) const
{
	return 0;
}

int Image32::Blur3X3(Image32& outputImage) const
{
	return 0;
}

int Image32::EdgeDetect3X3(Image32& outputImage) const
{
	return 0;
}
int Image32::ScaleNearest(const float& scaleFactor,Image32& outputImage) const
{
	outputImage.setSize(round(scaleFactor * width()), round(scaleFactor * height()));
	for (int u = 0; u < outputImage.width(); u++)
	for (int v = 0; v < outputImage.height(); v++)
	{
		outputImage.pixel(u, v) = NearestSample(u / scaleFactor, v / scaleFactor);
	}
	return 1;
}

int Image32::ScaleBilinear(const float& scaleFactor,Image32& outputImage) const
{
	outputImage.setSize(round(scaleFactor * width()), round(scaleFactor * height()));
	for (int u = 0; u < outputImage.width(); u++)
	for (int v = 0; v < outputImage.height(); v++)
	{
		outputImage.pixel(u, v) = BilinearSample(u / scaleFactor, v / scaleFactor);
	}
	return 1;
}

int Image32::ScaleGaussian(const float& scaleFactor,Image32& outputImage) const
{
	float radius = scaleFactor > 1 ? 1 : (1 / scaleFactor);
	float variance = pow(radius / 3, 2);

	outputImage.setSize(round(scaleFactor * width()), round(scaleFactor * height()));
	for (int u = 0; u < outputImage.width(); u++)
	for (int v = 0; v < outputImage.height(); v++)
	{
		outputImage.pixel(u, v) = GaussianSample(u / scaleFactor, v / scaleFactor, variance, radius);
	}
	return 1;
}

int Image32::RotateNearest(const float& angle,Image32& outputImage) const
{
	float centerX, centerY;
	SetRotatedSize(this, angle, outputImage, centerX, centerY);
	float radians = angle / 180 * PI;
	for (int u = 0; u < outputImage.width(); u++)
	for (int v = 0; v < outputImage.height(); v++)
	{
		float x = u - centerX;
		float y = v - centerY;

		float x_rot = x * cos(-radians) - y * sin(-radians);
		float y_rot = x * sin(-radians) + y * cos(-radians);

		outputImage.pixel(u, v) = NearestSample(x_rot, y_rot);
	}
	
	return 1;
}

int Image32::RotateBilinear(const float& angle,Image32& outputImage) const
{
	float centerX, centerY;
	SetRotatedSize(this, angle, outputImage, centerX, centerY);

	float radians = angle / 180 * PI;
	for (int u = 0; u < outputImage.width(); u++)
	for (int v = 0; v < outputImage.height(); v++)
	{
		float x = u - centerX;
		float y = v - centerY;

		float x_rot = x * cos(-radians) - y * sin(-radians);
		float y_rot = x * sin(-radians) + y * cos(-radians);

		outputImage.pixel(u, v) = BilinearSample(x_rot, y_rot);
	}

	return 1;
}
	
int Image32::RotateGaussian(const float& angle,Image32& outputImage) const
{
	float centerX, centerY;
	SetRotatedSize(this, angle, outputImage, centerX, centerY);

	float radians = angle / 180 * PI;
	for (int u = 0; u < outputImage.width(); u++)
	for (int v = 0; v < outputImage.height(); v++)
	{
		float x = u - centerX;
		float y = v - centerY;

		float x_rot = x * cos(-radians) - y * sin(-radians);
		float y_rot = x * sin(-radians) + y * cos(-radians);

		outputImage.pixel(u, v) = GaussianSample(x_rot, y_rot, 1.0 / 9, 1);
	}
	return 1;
}


int Image32::SetAlpha(const Image32& matte)
{
	for (int u = 0; u < width(); u++)
	for (int v = 0; v < height(); v++)
	{
		// Use maximum channel to specify alpha. In practice only one channel will
		// be used, so this will work regardless of whether it's r g or b.
		pixel(u, v).a = fmax(matte.pixel(u, v).r, fmax(matte.pixel(u, v).g, matte.pixel(u, v).b));
	}
	return 1;
}

int Image32::Composite(const Image32& overlay,Image32& outputImage) const
{
	outputImage.setSize(width(), height());
	for (int u = 0; u < width(); u++)
	for (int v = 0; v < height(); v++)
	{
		Pixel32& p = outputImage.pixel(u, v);
		const Pixel32& s = pixel(u, v);
		const Pixel32& o = overlay.pixel(u, v);

		float overlay_alpha = o.a / (float)255;
		float source_alpha = (1 - overlay_alpha) * (s.a / (float) 255);
		p.a = Clamp32((overlay_alpha + source_alpha) * 255);
		p.r = Clamp32(overlay_alpha * o.r + source_alpha * s.r );
		p.g = Clamp32(overlay_alpha * o.g + source_alpha * s.g);
		p.b = Clamp32(overlay_alpha * o.b + source_alpha * s.b);
	}

	return 1;
}

int Image32::CrossDissolve(const Image32& source,const Image32& destination,const float& blendWeight,Image32& ouputImage)
{
	ouputImage.setSize(source.width(), source.height());
	for (int u = 0; u < source.width(); u++)
	for (int v = 0; v < source.height(); v++)
	{
		Pixel32& p = ouputImage.pixel(u, v);
		const Pixel32& s = source.pixel(u, v);
		const Pixel32& d = destination.pixel(u, v);

		p.a = Clamp32(s.a * (1 - blendWeight) + d.a * blendWeight);
		p.r = Clamp32(s.r * (1 - blendWeight) + d.r * blendWeight);
		p.g = Clamp32(s.g * (1 - blendWeight) + d.g * blendWeight);
		p.b = Clamp32(s.b * (1 - blendWeight) + d.b * blendWeight);
	}
	return 1;
}
int Image32::Warp(const OrientedLineSegmentPairs& olsp,Image32& outputImage) const
{
	outputImage.setSize(width(), height());
	for (int u = 0; u < width(); u++)
	for (int v = 0; v < height(); v++)
	{
		float sourceX, sourceY;
		olsp.getSourcePosition(u, v, sourceX, sourceY);
		Pixel32 warped;
		
		if (sourceX < 0 || sourceX >= width() - 1 || sourceY < 0 || sourceY >= height() - 1)
		{
			// If sample is out of image bounds, use nearest pixel
			if (sourceX >= 0 && sourceX < width() - 1)
				warped = BilinearSample(sourceX, sourceY < 0 ? 0 : height() - 1);
			else if (sourceY >= 0 && sourceY < height() - 1)
				warped = BilinearSample(sourceX < 0 ? 0 : width() - 1, sourceY);
			else
				warped = BilinearSample(sourceX < 0 ? 0 : width() - 1, sourceY < 0 ? 0 : height() - 1);
		}
		else
		{
			warped = BilinearSample(sourceX, sourceY);
		}

		Pixel32& p = outputImage.pixel(u, v);
		p.a = warped.a;
		p.r = warped.r;
		p.g = warped.g;
		p.b = warped.b;
	}
	return 1;
}

int Image32::FunFilter(Image32& outputImage) const
{
	return 0;
}
int Image32::Crop(const int& x1,const int& y1,const int& x2,const int& y2,Image32& outputImage) const
{
	outputImage.setSize(labs(x1 - x2) + 1, labs(y1 - y2) + 1);
	for (int u = 0; u < outputImage.width(); u++)
	for (int v = 0; v < outputImage.height(); v++)
	{
		outputImage.pixel(u, v) = pixel(fmin(x1, x2) + u, fmin(y1, y2) + v);
	}
	return 1;
}


Pixel32 Image32::NearestSample(const float& x,const float& y) const
{
	Pixel32 pixel_out;
	const Pixel32& pixel_nearest = pixel(Mirror(round(x), width()), Mirror(round(y), height()));

	pixel_out.a = pixel_nearest.a;
	pixel_out.r = pixel_nearest.r;
	pixel_out.g = pixel_nearest.g;
	pixel_out.b = pixel_nearest.b;

	return pixel_out;
}
Pixel32 Image32::BilinearSample(const float& x,const float& y) const
{
	Pixel32 pixel_out;

	int x1 = floor(x);
	int x2 = x1 + 1;
	int y1 = floor(y);
	int y2 = y1 + 1;

	// Pixels outside the image are mirrored
	int x1_m = Mirror(x1, width());
	int x2_m = Mirror(x2, width());
	int y1_m = Mirror(y1, height());
	int y2_m = Mirror(y2, height());


	float dx = x - x1;

	Pixel32 a;
	a.a = pixel(x1_m, y1_m).a * (1 - dx) + pixel(x2_m, y1_m).a * dx;
	a.r = pixel(x1_m, y1_m).r * (1 - dx) + pixel(x2_m, y1_m).r * dx;
	a.g = pixel(x1_m, y1_m).g * (1 - dx) + pixel(x2_m, y1_m).g * dx;
	a.b = pixel(x1_m, y1_m).b * (1 - dx) + pixel(x2_m, y1_m).b * dx;

	Pixel32 b;
	b.a = pixel(x1_m, y2_m).a * (1 - dx) + pixel(x2_m, y2_m).a * dx;
	b.r = pixel(x1_m, y2_m).r * (1 - dx) + pixel(x2_m, y2_m).r * dx;
	b.g = pixel(x1_m, y2_m).g * (1 - dx) + pixel(x2_m, y2_m).g * dx;
	b.b = pixel(x1_m, y2_m).b * (1 - dx) + pixel(x2_m, y2_m).b * dx;

	float dy = y - y1;

	pixel_out.a = Clamp32(a.a * (1 - dy) + b.a * dy);
	pixel_out.r = Clamp32(a.r * (1 - dy) + b.r * dy);
	pixel_out.g = Clamp32(a.g * (1 - dy) + b.g * dy);
	pixel_out.b = Clamp32(a.b * (1 - dy) + b.g * dy);

	return pixel_out;
}
Pixel32 Image32::GaussianSample(const float& x,const float& y,const float& variance,const float& radius) const
{
	float weight_sum = 0;
	Pixel pixel_out;
	pixel_out.a = 0;
	pixel_out.b = 0;
	pixel_out.g = 0;
	pixel_out.r = 0;

	float window_top = ceil(y + radius);
	float window_bot = floor(y - radius);
	float window_right = ceil(x + radius);
	float window_left = floor(x - radius);

	for (int i = window_left; i <= window_right; i++)
	for (int j = window_bot; j <= window_top; j++)
	{
		float dist = Distance(i, j, x, y);
		if (dist > radius) continue;

		float weight = Gaussian(variance, dist);
		weight_sum += weight;

		const Pixel32& p = pixel(Mirror(i, width()), Mirror(j, height()));
		pixel_out.a += p.a * weight / 255;
		pixel_out.r += p.r * weight / 255;
		pixel_out.g += p.g * weight / 255;
		pixel_out.b += p.b * weight / 255;
	}

	pixel_out.a /= weight_sum;
	pixel_out.r /= weight_sum;
	pixel_out.g /= weight_sum;
	pixel_out.b /= weight_sum;

	return Pixel32(pixel_out);
}

double Gaussian(const float& variance, const double& x)
{
	return 1 / sqrt(2 * variance * PI)  *  exp(-pow(x, 2) / (2 * variance));
}

float Distance(const float& x1, const float& y1, const float& x2, const float& y2)
{
	return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
}

float Clamp32(const float& x)
{
	return fmax(0, fmin(x, 255));
}

float ClampF(const float& x)
{
	return fmax(0, fmin(x, 1));
}

float Mirror(const int& x, const int& cutoff)
{
	if (x >= 0)
	{
		if ((x / cutoff) % 2 == 1)
			return (cutoff - 1) - (x % cutoff);
		else
			return x % cutoff;
	}
	else
	{
		if ((-x / (cutoff + 1)) % 2 == 1)
			return cutoff + (x % -cutoff);
		else
			return -((x % cutoff) + 1);
	}
}

void SetRotatedSize(const Image32 * inputImage, const float& angle, Image32& outputImage, float& centerX, float& centerY)
{
	int width = inputImage->width();
	int height = inputImage->height();

	float radians = angle / 180 * PI;

	float corners_x[4];
	float corners_y[4];

	corners_x[0] = 0;
	corners_y[0] = 0;
	corners_x[1] = width - 1;
	corners_y[1] = 0;
	corners_x[2] = width - 1;
	corners_y[2] = height - 1;
	corners_x[3] = 0;
	corners_y[3] = height - 1;

	float rotated_x[4];
	float rotated_y[4];
	for (int i = 0; i < 4; i++) {
		rotated_x[i] = corners_x[i] * cos(radians) - corners_y[i] * sin(radians);
		rotated_y[i] = corners_x[i] * sin(radians) + corners_y[i] * cos(radians);
	}

	float minx = fmin(rotated_x[0], fmin(rotated_x[1], fmin(rotated_x[2], rotated_x[3])));
	float miny = fmin(rotated_y[0], fmin(rotated_y[1], fmin(rotated_y[2], rotated_y[3])));
	float maxx = fmax(rotated_x[0], fmax(rotated_x[1], fmax(rotated_x[2], rotated_x[3])));
	float maxy = fmax(rotated_y[0], fmax(rotated_y[1], fmax(rotated_y[2], rotated_y[3])));

	centerX = -minx;
	centerY = -miny;
	outputImage.setSize(maxx - minx + 1, maxy - miny + 1);
}
