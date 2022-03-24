//
// Created by indevn on 22.3.20.
//

#ifndef SMALLPT_FILM_HPP
#define SMALLPT_FILM_HPP


class Film {
public:
    Film (const Vector2 &resolution, const std::string &filename) :
            fullResolution(resolution),
            filename{filename},
            pixels{std::makeunique<Color[]>(Width() * Height())} {}
}

public:
    int Width() const { return (int) fullResolution.x; }
    int Height() const { return (int) fullResolution.y; }
    c


#endif //SMALLPT_FILM_HPP
