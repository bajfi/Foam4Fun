#pragma once

#include <vector.H>
#include <vector>

struct MapUtils
{
    static Foam::vector parseCells(Foam::Istream& is);
    static std::vector<std::vector<unsigned char>> mapMatrix(
      const std::vector<std::vector<unsigned char>>& matrix, const std::pair<int, int>& nxny
    ) noexcept;
    static std::vector<std::vector<unsigned char>> renderTextToMatrix(
      std::string_view text, std::string_view fontPath, float fontSize
    );
};
