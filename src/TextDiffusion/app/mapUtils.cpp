// Copyright (c) 2026 JackLee
//
// This software is released under the MIT License.
// https://opensource.org/licenses/MIT

#include "mapUtils.h"

// Force stb_truetype to use std:: versions of math functions
#define STBTT_sqrt(x) std::sqrt(x)
#define STBTT_pow(x, y) std::pow(x, y)
#define STBTT_cos(x) std::cos(x)
#define STBTT_acos(x) std::acos(x)

#define STB_TRUETYPE_IMPLEMENTATION
#include "stb_truetype.h"

#include <labelList.H>
#include <token.H>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

std::vector<std::vector<unsigned char>> MapUtils::renderTextToMatrix(
  std::string_view text, std::string_view fontPath, float fontSize
)
{
    // Load font file
    std::ifstream fontFile(fontPath.data(), std::ios::binary);
    if (!fontFile.is_open())
    {
        std::cerr << "Failed to open font: " << fontPath << "\n";
        return {};
    }
    std::vector<unsigned char> fontBuffer(
      (std::istreambuf_iterator<char>(fontFile)), std::istreambuf_iterator<char>()
    );
    if (fontBuffer.empty())
    {
        std::cerr << "Font buffer empty\n";
        return {};
    }

    stbtt_fontinfo font;
    if (!stbtt_InitFont(&font, fontBuffer.data(), 0))
    {
        std::cerr << "Failed to initialize font\n";
        return {};
    }

    // scaling and vertical metrics
    auto scale = stbtt_ScaleForPixelHeight(&font, fontSize);
    int ascent, descent, lineGap;
    stbtt_GetFontVMetrics(&font, &ascent, &descent, &lineGap);
    auto baseline = static_cast<int>(std::round(ascent * scale));
    auto height = static_cast<int>(std::round((ascent - descent) * scale));
    if (height <= 0)
        height = static_cast<int>(std::ceil(fontSize));

    // compute total width by summing advances (plus kerning)
    int width = 0;
    int text_size = text.size();
    for (auto i = 0; i < text_size; ++i)
    {
        int advance, lsb;
        int code = static_cast<unsigned char>(text[i]);
        stbtt_GetCodepointHMetrics(&font, code, &advance, &lsb);
        width += static_cast<int>(std::round(advance * scale));
        if (i + 1 < text_size)
        {
            int kern =
              stbtt_GetCodepointKernAdvance(&font, code, static_cast<unsigned char>(text[i + 1]));
            width += static_cast<int>(std::round(kern * scale));
        }
    }
    if (width <= 0)
        width = 1;

    // create pixel buffer (grayscale)
    std::vector<unsigned char> bitmap(width * height, 0);

    int xpos = 0;
    for (int i = 0; i < text_size; ++i)
    {
        int code = static_cast<unsigned char>(text[i]);
        int w, h, xoff, yoff;
        auto* glyph = stbtt_GetCodepointBitmap(&font, scale, scale, code, &w, &h, &xoff, &yoff);
        if (!glyph)
            continue;

        int destX = xpos + xoff;
        int destY = baseline + yoff;

        for (int yy = 0; yy < h; ++yy)
        {
            int dstRow = destY + yy;
            if (dstRow < 0 || dstRow >= height)
                continue;
            for (int xx = 0; xx < w; ++xx)
            {
                int dstCol = destX + xx;
                if (dstCol < 0 || dstCol >= width)
                    continue;
                bitmap[dstRow * width + dstCol] = glyph[yy * w + xx];
            }
        }

        stbtt_FreeBitmap(glyph, nullptr);

        int advance, lsb;
        stbtt_GetCodepointHMetrics(&font, code, &advance, &lsb);
        xpos += static_cast<int>(std::round(advance * scale));
        if (i + 1 < text_size)
        {
            int kern =
              stbtt_GetCodepointKernAdvance(&font, code, static_cast<unsigned char>(text[i + 1]));
            xpos += static_cast<int>(std::round(kern * scale));
        }
    }

    // reshape to 2D matrix
    std::vector<std::vector<unsigned char>> matrix(height, std::vector<unsigned char>(width, 0));
    for (auto y = 0; y < height; ++y)
    {
        for (auto x = 0; x < width; ++x)
        {
            matrix[y][x] = bitmap[y * width + x];
        }
    }
    return matrix;
}

Foam::vector MapUtils::parseCells(Foam::Istream& is)
{
    using namespace Foam;
    vector result;
    token t(is);
    word model(is);

    Info << "Parsing cells for model: " << model << endl;

    if (t.isPunctuation() && t.pToken() == token::BEGIN_LIST)
    {
        // Examine next token
        token t(is);

        // Optional zone name
        if (t.isWord())
        {
            // Examine next token
            is >> t;
        }
        is.putBack(t);

        if (t.isPunctuation())
        {
            // New-style: read a list of 3 values
            if (t.pToken() == token::BEGIN_LIST)
            {
                {
                    labelList n;
                    is >> n;
                }
                is >> result;
            }
            else
            {
                FatalIOErrorInFunction(is)
                  << "Incorrect token while reading n, expected '(', found " << t.info()
                  << exit(FatalIOError);
            }
        }
    }

    return result;
}

std::vector<std::vector<unsigned char>> MapUtils::mapMatrix(
  const std::vector<std::vector<unsigned char>>& matrix, const std::pair<int, int>& nxny
) noexcept
{
    auto [nx, ny] = nxny;
    if (nx <= 0 || ny <= 0)
    {
        return {};
    }

    const int srcHeight = matrix.size(), srcWidth = matrix.front().size();
    if (srcHeight == 0 || srcWidth == 0)
    {
        return std::vector<std::vector<unsigned char>>(ny, std::vector<unsigned char>(nx, 0));
    }

    if (srcWidth == nx && srcHeight == ny)
    {
        return matrix;
    }

    auto mappedMatrix =
      std::vector<std::vector<unsigned char>>(ny, std::vector<unsigned char>(nx, 0));

    const auto xScale = (nx == 1) ? 0.0 : static_cast<double>(srcWidth - 1) / (nx - 1);
    const auto yScale = (ny == 1) ? 0.0 : static_cast<double>(srcHeight - 1) / (ny - 1);

    for (int y = 0; y < ny; ++y)
    {
        const auto srcY = yScale * y;
        const auto y0 = static_cast<int>(std::floor(srcY));
        const auto y1 = std::min(y0 + 1, srcHeight - 1);
        const auto fy = srcY - y0;

        for (int x = 0; x < nx; ++x)
        {
            const auto srcX = xScale * (x);
            const auto x0 = static_cast<int>(std::floor(srcX));
            const auto x1 = std::min(x0 + 1, srcWidth - 1);
            const auto fx = srcX - x0;

            const auto v00 = static_cast<double>(matrix[y0][x0]);
            const auto v10 = static_cast<double>(matrix[y0][x1]);
            const auto v01 = static_cast<double>(matrix[y1][x0]);
            const auto v11 = static_cast<double>(matrix[y1][x1]);

            const auto v0 = v00 + (v10 - v00) * fx;
            const auto v1 = v01 + (v11 - v01) * fx;
            const auto v = v0 + (v1 - v0) * fy;

            const auto vi = std::lround(v);
            mappedMatrix[y][x] = static_cast<unsigned char>(std::clamp(vi, 0L, 255L));
        }
    }

    return mappedMatrix;
}
