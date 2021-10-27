#pragma once
#include "Vec3.h"

#include <iostream>
#include <fstream>
#include <vector>


class Image {
public:
	Image(int width, int height):
		p_width(width), p_height(height){
		p_pixels.resize(width * height);
	}
	inline int width() {return p_width;}
	inline int height() { return p_height; }
	void set_pixel(int x, int y, Vec3f color) {
		p_pixels[y * p_width + x] = color;
	}
	inline void fillBackground(const Vec3f& color = Vec3f(0.f, 0.f, 1.f)) {
		for (size_t y = 0; y < p_height; y++)
			for (size_t x = 0; x < p_width; x++) {
				Vec3f color0(0.f, 0.f, 1.f);
				Vec3f color1(1.f, 1.f, 1.f);
				float dy = static_cast<float> (y) / (p_height - 1);
				float alpha = (dy < 0.f) ? 0.f : (1.f < dy) ? 1.f : dy;
				set_pixel(x,y,mix(color1, color0, alpha));
			}
	}
	void savePPM(std::string path){
		std::ofstream out;
		out.open(path);
		out << "P3" << std::endl;
		out << p_width << " " << p_height << std::endl;
		out << 255 << std::endl;
		for (size_t y = 0; y < p_height; y++) {
			for (size_t x = 0; x < p_width; x++) {
				out << static_cast<unsigned int> (255.f * p_pixels[y * p_width + x][0]) << " "
					<< static_cast<unsigned int> (255.f * p_pixels[y * p_width + x][1]) << " "
					<< static_cast<unsigned int> (255.f * p_pixels[y * p_width + x][2]) << " ";
			}
			out << std::endl;
		}
		out.close();
	}


private:
	int p_width;
	int p_height;
	std::vector<Vec3f> p_pixels;


};