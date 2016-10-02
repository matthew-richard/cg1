#include "image.h"
#include <stdlib.h>
#include <math.h>

////////////////////////////
// Image processing stuff //
////////////////////////////
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
		outputImage.pixel(u, v).r = fmax(0, fmin(rand_r, 255));
		outputImage.pixel(u, v).r = fmax(0, fmin(rand_g, 255));
		outputImage.pixel(u, v).r = fmax(0, fmin(rand_b, 255));
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
			pixelOut.r = fmax(fmin(pixelIn.r * brightness, 255), 0);
			pixelOut.g = fmax(fmin(pixelIn.g * brightness, 255), 0);
			pixelOut.b = fmax(fmin(pixelIn.b * brightness, 255), 0);
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
			pixelOut.r = fmax(fmin(round(avgLuminance + contrast * (pixelIn.r - avgLuminance)), 255), 0);
			pixelOut.g = fmax(fmin(round(avgLuminance + contrast * (pixelIn.g - avgLuminance)), 255), 0);
			pixelOut.b = fmax(fmin(round(avgLuminance + contrast * (pixelIn.b - avgLuminance)), 255), 0);
		}
	}
	return 1;
}

int Image32::Saturate(const float& saturation,Image32& outputImage) const
{
	return 0;
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
	return 0;
}

int Image32::ScaleBilinear(const float& scaleFactor,Image32& outputImage) const
{
	return 0;
}

int Image32::ScaleGaussian(const float& scaleFactor,Image32& outputImage) const
{
	return 0;
}

int Image32::RotateNearest(const float& angle,Image32& outputImage) const
{
	return 0;
}

int Image32::RotateBilinear(const float& angle,Image32& outputImage) const
{
	return 0;
}
	
int Image32::RotateGaussian(const float& angle,Image32& outputImage) const
{
	return 0;
}


int Image32::SetAlpha(const Image32& matte)
{
	setSize(matte.width(), matte.height());

	for (int u = 0; u < width(); u++)
	for (int v = 0; v < height(); v++)
	{
		// Use maximum channel to specifiy alpha. In practice only one channel will
		// be used, so this will work regardless of whether it's r g or b.
		pixel(u, v).a = fmax(matte.pixel(u, v).r, fmax(matte.pixel(u, v).g, matte.pixel(u, v).b));
	}
	return 1;
}

int Image32::Composite(const Image32& overlay,Image32& outputImage) const
{
	return 0;
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

		p.a = s.a * (1 - blendWeight) + d.a * blendWeight;
		p.r = s.r * (1 - blendWeight) + d.r * blendWeight;
		p.g = s.g * (1 - blendWeight) + d.g * blendWeight;
		p.b = s.b * (1 - blendWeight) + d.b * blendWeight;
	}
	return 1;
}
int Image32::Warp(const OrientedLineSegmentPairs& olsp,Image32& outputImage) const
{
	// warp constants
	float a = 0.000001;
	float b = 0.5;
	float c = 0.9;

	for (int u = 0; u < width(); u++)
	for (int v = 0; v < height(); v++)
	{
		float sourceX, sourceY;
		olsp.getSourcePosition(u, v, sourceX, sourceY);
		Pixel32 warped = BilinearSample(sourceX, sourceY);

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
	return 0;
}


Pixel32 Image32::NearestSample(const float& x,const float& y) const
{
	Pixel32 pixel_out;

	if (x < -0.5 || x >= width() - 1 + 0.5 || y < -0.5 || y >= height() - 1 + 0.5)
	{
		pixel_out.a = 1;
		pixel_out.r = 0;
		pixel_out.g = 0;
		pixel_out.b = 0;
		return pixel_out;
	}

	const Pixel32& pixel_nearest = pixel(round(x), round(y));

	pixel_out.a = pixel_nearest.a;
	pixel_out.r = pixel_nearest.r;
	pixel_out.g = pixel_nearest.g;
	pixel_out.b = pixel_nearest.b;

	return pixel_out;
}
Pixel32 Image32::BilinearSample(const float& x,const float& y) const
{
	Pixel32 pixel_out;
	if (x < 0 || x >= width() - 1 || y < 0 || y >= height() - 1)
	{
		pixel_out.a = 1;
		pixel_out.r = 0;
		pixel_out.g = 0;
		pixel_out.b = 0;
		return pixel_out;
	}

	int x1 = floor(x);
	int x2 = x1 + 1;
	int y1 = floor(y);
	int y2 = y1 + 1;

	float dx = x - x1;

	Pixel32 a;
	a.a = pixel(x1, y1).a * (1 - dx) + pixel(x2, y1).a * dx;
	a.r = pixel(x1, y1).r * (1 - dx) + pixel(x2, y1).r * dx;
	a.g = pixel(x1, y1).g * (1 - dx) + pixel(x2, y1).g * dx;
	a.b = pixel(x1, y1).b * (1 - dx) + pixel(x2, y1).b * dx;

	Pixel32 b;
	b.a = pixel(x1, y2).a * (1 - dx) + pixel(x2, y2).a * dx;
	b.r = pixel(x1, y2).r * (1 - dx) + pixel(x2, y2).r * dx;
	b.g = pixel(x1, y2).g * (1 - dx) + pixel(x2, y2).g * dx;
	b.b = pixel(x1, y2).b * (1 - dx) + pixel(x2, y2).b * dx;

	float dy = y - y1;

	pixel_out.a = a.a * (1 - dy) + b.a * dy;
	pixel_out.r = a.r * (1 - dy) + b.r * dy;
	pixel_out.g = a.g * (1 - dy) + b.g * dy;
	pixel_out.b = a.b * (1 - dy) + b.g * dy;

	return pixel_out;
}
Pixel32 Image32::GaussianSample(const float& x,const float& y,const float& variance,const float& radius) const
{
	float weight_sum = 0;
	Pixel32 pixel_out;
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
		if (i < 0 || i >= width() - 1 || j < 0 || j >= height() - 1) continue;

		float dist = Distance(i, j, x, y);
		if (dist > radius) continue;

		float weight = Gaussian(variance, dist);
		weight_sum += weight;

		pixel_out.a += pixel(i, j).a * weight;
		pixel_out.r += pixel(i, j).r * weight;
		pixel_out.g += pixel(i, j).g * weight;
		pixel_out.b += pixel(i, j).b * weight;
	}

	pixel_out.a /= weight_sum;
	return pixel_out;
}

double Gaussian(const float& variance, const double& x)
{
	return 1 / sqrt(2 * variance * PI)  *  exp(-pow(x, 2) / (2 * variance));
}

double Distance(const float& x1, const float&& y1, const float& x2, const float& y2)
{
	return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
}
