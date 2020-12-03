#include <list>
#include <SDL.h>
#include <chrono>
#include <ctime>  
#include <iostream>
#include <limits>
#include <random>
#include <vector>
#include <algorithm>
#include <math.h>
static std::mt19937 rnd;
#define random_color() std::uniform_int_distribution<>(0x000001, 0xFFFFFF)(rnd)
#define random_coord(max) std::uniform_int_distribution<>(0, max - 1)(rnd)
const int FPS = 30;
int loopcounter = 0;
const int FRAMETIME = 1000 / FPS;
unsigned int screen_w = 200;
unsigned int screen_h = 200;
void greetings()
{
    std::cout << "press one of buttons below to do things" << std::endl;
    std::cout << "LMB on screen: set seed in place of cursor" << std::endl;

    std::cout << "1: homogeneous seed distribution" << std::endl;
    std::cout << "2: random seed distribution" << std::endl;
    std::cout << "3: random seed distribution with minimum range between" << std::endl;
    std::cout << "4: change grid size" << std::endl;
    std::cout << "5: switch periodical/nonperiodical" << std::endl;
    std::cout << "6: change nucleation mode" << std::endl;
    std::cout << "7: toggle monte carlo" << std::endl;

}
bool setPixel(SDL_Surface* surface, int x, int y, Uint32 pixel)
{
    Uint32* target_pixel = (Uint32*)((Uint8*)surface->pixels + y * surface->pitch + x * sizeof * target_pixel);
    if (target_pixel == 0) return 1;
    *target_pixel = pixel;
    return 0;
}
inline Uint32 probePixel(SDL_Surface* surface, int x, int y, bool period)
{
    //Uint32 target_pixel;
    //if (period)
    //{
    //    target_pixel = *(Uint32*)((Uint8*)surface->pixels + (y % screen_h) * surface->pitch + (x % screen_w) * sizeof * &target_pixel);
    //}
    //else
    //{
    //    if (x < 0 || y < 0 || x >= screen_w || y >= screen_h) return 0;
    //    target_pixel = *(Uint32*)((Uint8*)surface->pixels + y * surface->pitch + x * sizeof * &target_pixel);
    //}
    if (period)
    {

    }
    else
    {
        if (x < 0 || y < 0 || x >= (int)screen_w || y >= (int)screen_h) return 0;
        //target_pixel = ((Uint32*)surface->pixels)[y * surface->w + x];
        int bpp = surface->format->BytesPerPixel;
        /* Here p is the address to the pixel we want to retrieve */
        Uint8* p = (Uint8*)surface->pixels + y * surface->pitch + x * bpp;
        return *(Uint32*)p;
    }
    //return target_pixel;
}
struct Point {
    int x;
    int y;
    bool operator==(const Point p)
    {
        if (this->x == p.x && this->y == p.y) return true;
        else return false;
    }
    Point(int _x, int _y) { x = _x; y = _y; }
    bool isInRadius(const Point* p, int r)
    {
        double distance = sqrt(pow(this->x - p->x, 2) + pow(this->y - p->y, 2));
        if (distance > r) return false;
        else return true;
    }
};
std::vector<int> color_register;
std::vector<Point> borders;
void setPixelRandomColor(SDL_Renderer* renderer, SDL_Surface* surface, int x, int y)
{
    if (probePixel(surface, x, y, 0) != 0) return;
    int pixel_color;
    bool loop;
    do {
        pixel_color = random_color();
        loop = !all_of(color_register.begin(), color_register.end(), [](int i) { return i != 0; });
        setPixel(surface, x, y, pixel_color);

    } while (loop);

    setPixel(surface, x, y, pixel_color);
    color_register.push_back(pixel_color);
    auto texture = SDL_CreateTextureFromSurface(renderer, surface);
    SDL_RenderCopy(renderer, texture, NULL, NULL);
    SDL_RenderPresent(renderer);
    SDL_DestroyTexture(texture);
}
void setPixelRandomColor2(SDL_Surface* surface, int x, int y)
{
    if (probePixel(surface, x, y, 0) != 0) return;
    int pixel_color;
    bool loop;
    do {
        pixel_color = random_color();
        loop = !all_of(color_register.begin(), color_register.end(), [](int i) { return i != 0; });
        setPixel(surface, x, y, pixel_color);

    } while (loop);

    setPixel(surface, x, y, pixel_color);
    color_register.push_back(pixel_color);
}

enum nucleationModes {
    VONNEUMAN, MOORE, PENTALEFT, PENTARIGHT, PENTAUP, PENTADOWN, PENTARANDOM, HEXALEFT, HEXARIGHT, HEXARANDOM, MONTECARLO
};
int nucleation_mode = HEXARANDOM;
int recent_nucleation_mode = HEXARANDOM;
bool periodicalness;
void vonNeumann(SDL_Surface* s, int x, int y, std::vector<std::pair<int, int>>* neighbors, bool period)
{
    int pxc[4];
    pxc[0] = probePixel(s, x + 1, y, period);
    pxc[1] = probePixel(s, x - 1, y, period);
    pxc[2] = probePixel(s, x, y + 1, period);
    pxc[3] = probePixel(s, x, y - 1, period);
    for (int ii = 0; ii < 4; ii++)
    {
        if (pxc[ii] == 0) continue;
        std::pair<int, int> to_push;
        to_push.first = pxc[ii];
        to_push.second++;
        neighbors->push_back(to_push);
    }
}
void moore(SDL_Surface* s, int x, int y, std::vector<std::pair<int, int>>* neighbors, bool period)
{
    for (int i = y - 1; i <= y + 1; i++)
        for (int j = x - 1; j <= x + 1; j++)
        {
            if (i == y && j == x) continue;
            int pxc = probePixel(s, j, i, period);
            if (pxc == 0) continue;
            std::pair<int, int> to_push;
            to_push.first = pxc;
            to_push.second++;
            neighbors->push_back(to_push);
        }
}
void pentaLeft(SDL_Surface* s, int x, int y, std::vector<std::pair<int, int>>* neighbors, bool period)
{
    int pxc[5];
    pxc[0] = probePixel(s, x, y + 1, period);
    pxc[1] = probePixel(s, x, y - 1, period);
    pxc[2] = probePixel(s, x - 1, y, period);
    pxc[3] = probePixel(s, x - 1, y + 1, period);
    pxc[4] = probePixel(s, x - 1, y - 1, period);
    for (int ii = 0; ii < 5; ii++)
    {
        if (pxc[ii] == 0) continue;
        std::pair<int, int> to_push;
        to_push.first = pxc[ii];
        to_push.second++;
        neighbors->push_back(to_push);
    }
}
void pentaRight(SDL_Surface* s, int x, int y, std::vector<std::pair<int, int>>* neighbors, bool period)
{
    int pxc[5];
    pxc[0] = probePixel(s, x, y + 1, period);
    pxc[1] = probePixel(s, x, y - 1, period);
    pxc[2] = probePixel(s, x + 1, y, period);
    pxc[3] = probePixel(s, x + 1, y + 1, period);
    pxc[4] = probePixel(s, x + 1, y - 1, period);
    for (int ii = 0; ii < 5; ii++)
    {
        if (pxc[ii] == 0) continue;
        std::pair<int, int> to_push;
        to_push.first = pxc[ii];
        to_push.second++;
        neighbors->push_back(to_push);
    }
}
void pentaUp(SDL_Surface* s, int x, int y, std::vector<std::pair<int, int>>* neighbors, bool period)
{
    int pxc[5];
    pxc[0] = probePixel(s, x, y + 1, period);
    pxc[1] = probePixel(s, x + 1, y + 1, period);
    pxc[2] = probePixel(s, x - 1, y + 1, period);
    pxc[3] = probePixel(s, x - 1, y, period);
    pxc[4] = probePixel(s, x + 1, y, period);
    for (int ii = 0; ii < 5; ii++)
    {
        if (pxc[ii] == 0) continue;
        std::pair<int, int> to_push;
        to_push.first = pxc[ii];
        to_push.second++;
        neighbors->push_back(to_push);
    }
}
void pentaDown(SDL_Surface* s, int x, int y, std::vector<std::pair<int, int>>* neighbors, bool period)
{
    int pxc[5];
    pxc[0] = probePixel(s, x, y - 1, period);
    pxc[1] = probePixel(s, x + 1, y - 1, period);
    pxc[2] = probePixel(s, x - 1, y - 1, period);
    pxc[3] = probePixel(s, x - 1, y, period);
    pxc[4] = probePixel(s, x + 1, y, period);
    for (int ii = 0; ii < 5; ii++)
    {
        if (pxc[ii] == 0) continue;
        std::pair<int, int> to_push;
        to_push.first = pxc[ii];
        to_push.second++;
        neighbors->push_back(to_push);
    }
}
void hexaLeft(SDL_Surface* s, int x, int y, std::vector<std::pair<int, int>>* neighbors, bool period)
{
    for (int i = y - 1; i <= y + 1; i++)
        for (int j = x - 1; j <= x + 1; j++)
        {
            if (i == y + 1 && j == x - 1 || i == y - 1 && j == x + 1) continue;
            int pxc = probePixel(s, j, i, period);
            if (pxc == 0) continue;
            std::pair<int, int> to_push;
            to_push.first = pxc;
            to_push.second++;
            neighbors->push_back(to_push);
        }
}
void hexaRight(SDL_Surface* s, int x, int y, std::vector<std::pair<int, int>>* neighbors, bool period)
{
    for (int i = y - 1; i <= y + 1; i++)
        for (int j = x - 1; j <= x + 1; j++)
        {
            if (i == y - 1 && j == x - 1 || i == y + 1 && j == x + 1) continue;
            int pxc = probePixel(s, j, i, period);
            if (pxc == 0) continue;
            std::pair<int, int> to_push;
            to_push.first = pxc;
            to_push.second++;
            neighbors->push_back(to_push);
        }
}
Uint32 neighborsAndProbe(SDL_Surface* s, int x, int y, std::vector<std::pair<int, int>>* neighbors, bool period)
{
    Uint32 pixel_color = probePixel(s, x, y, period);
    for (int i = y - 1; i <= y + 1; i++)
        for (int j = x - 1; j <= x + 1; j++)
        {
            if (i == y && j == x) continue;
            int pxc = probePixel(s, j, i, period);
            if (pxc == 0 /*|| pxc == pixel_color*/) continue;
            std::pair<int, int> to_push;
            to_push.first = pxc;
            to_push.second = 1;
            neighbors->push_back(to_push);
        }
    return pixel_color;
}
inline void searchForBorders(SDL_Surface* s)
{
    for (int y = 0; y < s->h; y++)
        for (int x = 0; x < s->w; x++)
        {
            int pxc = probePixel(s, x, y, 0);
            if (pxc == 0) continue;
            std::vector < std::pair<int, int>> n;
            moore(s, x, y, &n, 0);
            if (any_of(n.begin(), n.end(), [pxc](std::pair<int, int> p) {if (p.first != pxc) return true; else return false; }))
            {
                Point p(x, y);
                borders.push_back(p);
            }
        }
}
//inline void monteCarlo(SDL_Surface* s, bool period)
//{
//    //if (borders.empty())
//    //{
//    std::cout << loopcounter++ << std::endl;
//    searchForBorders(s);
//    //}
//    for (int i = 0; i < s->h; i++)
//        for (int j = 0; j < s->h; j++)
//        {
//            //int n = random_coord(borders.size());
//            //Point pointeuuu = borders[n];
//            //int ptc = probePixel(s, borders[n].x, borders[n].y, 0);
//            //int ptc = probePixel(s, j, i, 0);
//            //int ptcr;
//            std::vector<std::pair<int, int>> neighbors;
//            int ptc = neighborsAndProbe(s, j, i, &neighbors, 0);
//            int energy1 = 0;
//            int energy2 = 0;
//            for (std::pair<int, int> pair : neighbors)
//            {
//                if (pair.first != ptc) energy1 += pair.second;
//
//            }
//
//            int ptcr = neighbors[random_coord(neighbors.size())].first;
//            for (std::pair<int, int> pair : neighbors)
//            {
//                if (pair.first != ptcr) energy2 += pair.second;
//
//            }
//            if (energy2 < energy1)
//            {
//                setPixel(s, j, i, ptcr);
//            }
//            //borders.erase(borders.begin() + n);
//            //borders.shrink_to_fit();
//            //for(int i = 0; i < borders.size())
//        }
//}
inline void monteCarlo(SDL_Surface* s, bool period)
{
    for (int i = 0; i < s->h; i++)
    {
        for (int j = 0; j < s->w; j++)
        {
            std::vector<std::pair<int, int>> neighbors;
            Uint32 pixel_color = neighborsAndProbe(s, j, i, &neighbors, periodicalness);
            if (neighbors.empty())
                continue;
            int original_energy = 0;
            for (auto p : neighbors)
            {
                if (p.first != pixel_color)
                    original_energy++;
            }
            Uint32 new_pixel_color = neighbors[random_coord(neighbors.size())].first;
            int new_energy = 0;
            for (auto p : neighbors)
            {
                if (p.first != new_pixel_color)
                    new_energy++;
            }
            if (new_energy <= original_energy)
            {
                setPixel(s, j, i, new_pixel_color);
            }
        }
    }
}
bool isBorder(SDL_Surface* s, int x, int y, bool period)
{
    int pxlc = probePixel(s, x, y, period);
    for (int i = y - 1; i <= y + 1; i++)
        for (int j = x - 1; j <= x + 1; j++)
        {
            if (i == y && j == x) continue;
            int pxc = probePixel(s, j, i, period);
            if (pxc == 0) continue;
            if (pxc != pxlc) return true;
        }
    return false;
}
//inline bool recrystalise(SDL_Window *window, SDL_Renderer *renderer)
//{
//    //all the declarations and initialisations
//    auto windowSurface = SDL_GetWindowSurface(window);
//    static bool recrystalisation_finished = false;
//    static double THE_MYSTERY_VALUE = 0;
//    static SDL_Surface* recrystalisation_reference_surface = SDL_CreateRGBSurface(NULL, windowSurface->w, windowSurface->h, 32, 0, 0, 0, 0);
//    static SDL_Surface* crystalised_structure = SDL_CreateRGBSurface(NULL, windowSurface->w, windowSurface->h, 32, 0, 0, 0, 0);
//    static SDL_Surface* crystalised_structure_copy = SDL_CreateRGBSurface(NULL, windowSurface->w, windowSurface->h, 32, 0, 0, 0, 0);
//    static SDL_Surface* prenucleidation = SDL_CreateRGBSurface(NULL, windowSurface->w, windowSurface->h, 32, 0, 0, 0, 0);
//    SDL_SetColorKey(crystalised_structure, SDL_TRUE, 0x00000000);
//    SDL_SetColorKey(prenucleidation, SDL_TRUE, 0x00000000);
//    SDL_BlitSurface(crystalised_structure, NULL, crystalised_structure_copy, NULL);
//    static double time = 0;
//    double LAST_MYSTERY_VALUE = THE_MYSTERY_VALUE;
//    THE_MYSTERY_VALUE = 86710969050178.5 / 9.41268203527779 + (1 - 86710969050178.5 / 9.41268203527779) * exp(-9.41268203527779 * time);//ro
//    double THE_MYSTERY_VALUE_CRITICAL_VALUE_OVERLOAD = 1.9 * 0.000000000257 * 8 *pow(10, 10) * sqrt(time);//sigma
//
//    static bool dislocation_table_initialized = false;
//    static double **dislocation_table; //[x][y]
//    static double global_dislocation_value_border = 0;
//    static double global_dislocation_value_interior = 0;
//
//    double THE_MYSTERY_VALUE_ATTRIBUTED_TO_THE_AREA = (THE_MYSTERY_VALUE - LAST_MYSTERY_VALUE);
//    global_dislocation_value_border += 0.7 * THE_MYSTERY_VALUE_ATTRIBUTED_TO_THE_AREA;
//    global_dislocation_value_interior += 0.3 * THE_MYSTERY_VALUE_ATTRIBUTED_TO_THE_AREA;
//    double asserted_dislocation_total = 0;
//    double asserted_dislocation_border = 0;
//    double asserted_dislocation_interior = 0;
//    double THE_MYSTERY_VALUE_CHUNK = THE_MYSTERY_VALUE_ATTRIBUTED_TO_THE_AREA / (int)((int)windowSurface->w * (int)windowSurface->h);
//    double asserted_dislocation_border_max = 0.8 * THE_MYSTERY_VALUE_ATTRIBUTED_TO_THE_AREA;
//    double asserted_dislocation_interior_max = 0.2 * THE_MYSTERY_VALUE_ATTRIBUTED_TO_THE_AREA;
//
//    // przedrozrost
//    for (int y = 0; y < windowSurface->h; y++)
//        for (int x = 0; x < windowSurface->w; x++)
//        {
//            std::vector<std::pair<int, int>> neighbors;
//            int pxc = neighborsAndProbe(crystalised_structure, x, y, &neighbors, periodicalness);
//            if (pxc != 0) continue;
//            if (dislocation_table[x][y + 1] < dislocation_table[x][y] &&
//                dislocation_table[x][y - 1] < dislocation_table[x][y] &&
//                dislocation_table[x + 1][y] < dislocation_table[x][y] &&
//                dislocation_table[x - 1][y] < dislocation_table[x][y])
//                setPixel(prenucleidation, x, y, );
//        }
//    // rozrost
//    for (int i = 0; i < windowSurface->w; i++)
//        for (int j = 0; j < windowSurface->h; j++)
//        {
//            Uint32 pxlc = probePixel(crystalised_structure, j, i, periodicalness);
//            if (pxlc != 0) continue;
//
//            std::vector<std::pair<int, int>> neighbors;
//            neighborsAndProbe(crystalised_structure, j, i, &neighbors, periodicalness);
//            if (neighbors.size() > 0)
//            {
//                setPixel(crystalised_structure_copy, j, i, neighbors[random_coord(neighbors.size())].first);
//            }
//        }
//    // rozrost koniec
//    SDL_BlitSurface(crystalised_structure_copy, NULL, crystalised_structure, NULL);
//    auto t = SDL_CreateTextureFromSurface(renderer, crystalised_structure);
//    SDL_RenderCopy(renderer, t, NULL, NULL);
//
//    SDL_RenderPresent(renderer);
//
//    if (!dislocation_table_initialized)
//    {
//        recrystalisation_reference_surface = SDL_CreateRGBSurface(NULL, windowSurface->w, windowSurface->h, 32, NULL, NULL, NULL, NULL);
//        SDL_BlitSurface(windowSurface, NULL, recrystalisation_reference_surface, NULL);
//        dislocation_table = new double* [windowSurface->w];
//        for (int i = 0; i < windowSurface->w; i++)
//        {
//            dislocation_table[i] = new double[windowSurface->h];
//            for (int j = 0; j < windowSurface->h; j++)
//                dislocation_table[i][j] = 0;
//        }
//        dislocation_table_initialized = true;
//    }
//
//    //robienie rzeczy
//
//    // Rozrost
//
//    // Update dislocation
//    while (asserted_dislocation_interior < asserted_dislocation_interior_max && asserted_dislocation_border < asserted_dislocation_border_max)
//    {
//        int x = random_coord(windowSurface->w);
//        int y = random_coord(windowSurface->w);
//        if (probePixel(recrystalisation_reference_surface, x, y, periodicalness) == 0) continue;
//        bool border = isBorder(recrystalisation_reference_surface, x, y, periodicalness);
//        if (!border && asserted_dislocation_interior >= asserted_dislocation_interior_max || border && asserted_dislocation_border >= asserted_dislocation_border_max) continue;
//        dislocation_table[x][y] += THE_MYSTERY_VALUE_CHUNK;
//        
//        if (border) asserted_dislocation_border += THE_MYSTERY_VALUE_CHUNK;
//        else asserted_dislocation_interior += THE_MYSTERY_VALUE_CHUNK;
//        asserted_dislocation_total += THE_MYSTERY_VALUE_CHUNK;
//        if (border)
//        {
//            if(dislocation_table[x][y] + global_dislocation_value_border >= THE_MYSTERY_VALUE_ATTRIBUTED_TO_THE_AREA)
//                setPixelRandomColor(renderer, crystalised_structure, x, y);
//        }
//        else
//        {
//            if (dislocation_table[x][y] + global_dislocation_value_interior >= THE_MYSTERY_VALUE_ATTRIBUTED_TO_THE_AREA)
//                setPixelRandomColor(renderer, crystalised_structure, x, y);
//        }
//    }
//    //auto recrystalisation_reference_surface_texture = SDL_CreateTextureFromSurface(renderer, recrystalisation_reference_surface);
//    //auto crystalised_structure_texture = SDL_CreateTextureFromSurface(renderer, crystalised_structure);
//    //SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);
//    //SDL_SetTextureBlendMode(recrystalisation_reference_surface_texture, SDL_BLENDMODE_BLEND);
//
//    //SDL_Surface* output_surface = SDL_CreateRGBSurface(NULL, windowSurface->w, windowSurface->h, 0, 0, 0, 0, 0);
//    //SDL_BlitSurface(recrystalisation_reference_surface, NULL, output_surface, NULL);
//    //SDL_BlitSurface(crystalised_structure, NULL, output_surface, NULL);
//
//    auto output_texture = SDL_CreateTextureFromSurface(renderer, crystalised_structure);
//
//    SDL_RenderCopy(renderer, output_texture, NULL, NULL);
//
//    SDL_RenderPresent(renderer);
//
//    //SDL_FreeSurface(output_surface);
//    SDL_DestroyTexture(output_texture);
//
//    //koniec robienia rzeczy
//
//    if (recrystalisation_finished)
//    {
//        SDL_FreeSurface(recrystalisation_reference_surface);
//        SDL_FreeSurface(crystalised_structure);
//    }
//
//    std::cout << time << std::endl;
//    time = time + 0.001;
//
//    return recrystalisation_finished;
//}
inline bool recrystalise(SDL_Window* window, SDL_Renderer* renderer)
{
    static double time = 0;
    static bool initialised = false;
    static bool recrystalisation_finished = false;
    static double **dislocation_table;
    double dislocation;
    double dislocation_per_cell;
    double dislocation_new;
    static double dislocation_old;
    double dislocation_limit;
    double dislocation_package;
    double dislocation_leftover;
    double dislocation_chunk;
    static double dislocation_global_interior;
    static double dislocation_global_border;
    static SDL_Surface* recrystalisation_structure;
    SDL_Surface* recrystalisation_structure_new;
    static SDL_Surface* recrystalisation_base;
    static int w;
    static int h;

    // Initialisation
    if (!initialised)
    {
        dislocation_global_interior = 0;
        dislocation_global_border = 0;
        SDL_Surface* windowSurface = SDL_GetWindowSurface(window);
        w = windowSurface->w;
        h = windowSurface->h;
        recrystalisation_structure = SDL_CreateRGBSurface(NULL, w, h, 32, 0, 0, 0, 0);
        SDL_SetColorKey(recrystalisation_structure, SDL_TRUE, 0x00000000);
        //recrystalisation_base = SDL_CreateRGBSurface(NULL, w, h, 32, 0, 0, 0, 0);
        recrystalisation_base = SDL_DuplicateSurface(windowSurface);

        dislocation_table = new double*[w];
        for (int x = 0; x < w; x++)
        {
            dislocation_table[x] = new double[h];
            for (int y = 0; y < h; y++)
                dislocation_table[x][y] = 0;
        }

        initialised = true;
    }
    recrystalisation_structure_new = SDL_CreateRGBSurface(NULL, w, h, 32, 0, 0, 0, 0);
    SDL_SetColorKey(recrystalisation_structure_new, SDL_TRUE, 0x00000000);

    // Update variables
    dislocation_new = 86710969050178.5 / 9.41268203527779 + (1 - 86710969050178.5 / 9.41268203527779) * exp(-9.41268203527779 * time);
    dislocation = dislocation_new - dislocation_old;
    //tension??? = 1.9 * 86710969050178.5 * 9.41268203527779 * sqrt(dislocation);
    dislocation_per_cell = dislocation / (static_cast<int>(w) * static_cast<int>(h));

    // Nucleate
    for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++)
        {
            Uint32 pxlc = probePixel(recrystalisation_structure, x, y, periodicalness);
            if (pxlc != 0) continue;

            std::vector<std::pair<int, int>> neighbors;
            neighborsAndProbe(recrystalisation_structure, x, y, &neighbors, periodicalness);
            if (neighbors.size() > 0)
            {
                setPixel(recrystalisation_structure_new, x, y, neighbors[random_coord(neighbors.size())].first);
            }
        }

    // Distribute dislocation
    dislocation_leftover = 0;
    for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++)
        {
            if(probePixel(recrystalisation_structure, x, y, periodicalness) != 0) continue;
            if (isBorder(recrystalisation_base, x, y, periodicalness))
            {
                dislocation_table[x][y] += 0.7 * dislocation_per_cell;
                dislocation_leftover += 0.3 * dislocation_per_cell;
            }
            else
            {
                dislocation_table[x][y] += 0.3 * dislocation_per_cell;
                dislocation_leftover += 0.7 * dislocation_per_cell;
            }

        }

    dislocation_chunk = dislocation_leftover / 100;//static_cast<int>(sqrt(w * h));
    while (dislocation_leftover > 0)
    {
        int x = random_coord(w);
        int y = random_coord(h);

        if (probePixel(recrystalisation_structure, x, y, periodicalness) != 0) continue;

        if (isBorder(recrystalisation_base, x, y, periodicalness))
        {
            dislocation_table[x][y] += 0.8 * dislocation_chunk;
            dislocation_leftover -= 0.8 * dislocation_chunk;
        }
        else
        {
            dislocation_table[x][y] += 0.2 * dislocation_chunk;
            dislocation_leftover -= 0.2 * dislocation_chunk;
        }
    }

    // checking for new seeds + growth
    for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++)
        {
            if (probePixel(recrystalisation_structure, x, y, periodicalness) != 0) continue;

            //if (y + 1 > h ||
            //    y - 1 < 0 ||
            //    x + 1 > w ||
            //    x - 1 < 0) continue;
            //if ((y + 1 > h) ? true : (dislocation_table[x][y + 1] < dislocation_table[x][y]) &&
            //    (y - 1 < 0) ? true : (dislocation_table[x][y - 1] < dislocation_table[x][y]) &&
            //    (x + 1 > w) ? true : (dislocation_table[x + 1][y] < dislocation_table[x][y]) &&
            //    (x - 1 < 0) ? true : (dislocation_table[x - 1][y] < dislocation_table[x][y]))
            //    setPixelRandomColor2(recrystalisation_structure_new, x, y);
            //else
            if (dislocation_table[x][y] > (dislocation))
            {
                std::cout << "new seed at (" << x << "," << y << ")" << std::endl;
                setPixelRandomColor2(recrystalisation_structure_new, x, y);
            }
        }

    SDL_BlitSurface(recrystalisation_structure_new, NULL, recrystalisation_structure, NULL);
    SDL_Surface* windowSurface = SDL_GetWindowSurface(window);
    SDL_BlitSurface(recrystalisation_structure, NULL, windowSurface, NULL);
    SDL_Texture* t = SDL_CreateTextureFromSurface(renderer, windowSurface);
    SDL_RenderCopy(renderer, t, NULL, NULL);
    SDL_RenderPresent(renderer);
    time = time + 0.001;
    //std::cout << time << " " << dislocation / (**dislocation_table) << std::endl;
    SDL_FreeSurface(recrystalisation_structure_new);
    SDL_DestroyTexture(t);

    return recrystalisation_finished;
}

inline void mainLoop()
{
    greetings();
    bool quit = false;
    bool pause = true;
    bool montecarlo = false;
    bool recrystalisation = false;

    SDL_Event event;
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* window = SDL_CreateWindow("pixels on screen", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, screen_w, screen_h, SDL_WINDOW_RESIZABLE);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, 0);
    auto start = std::chrono::system_clock::now();
    int the_time = 1;

    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    auto windowSurface = SDL_GetWindowSurface(window);
    auto texture = SDL_CreateTextureFromSurface(renderer, windowSurface);
    SDL_RenderCopy(renderer, texture, NULL, NULL);
    bool lock_pause = false;
    bool lock_input_1 = false;
    bool lock_input_2 = false;
    bool lock_input_3 = false;
    SDL_RenderPresent(renderer);

    while (!quit)
    {
        bool click = false;
        SDL_PollEvent(&event);
        if (event.type != SDL_MOUSEBUTTONUP) click = false;
        //if (event.key.repeat == true) continue;
        switch (event.type)
        {
        case SDL_KEYDOWN:
            switch (event.key.keysym.sym)
            {
            case SDLK_SPACE:
                if (!pause && !lock_pause)
                {
                    // SPC D, P0, PP0
                    pause = true;
                    lock_pause = true;
                    std::cout << "paused and locked pause" << std::endl;
                }
                else if (pause && !lock_pause)
                {
                    pause = false;
                    lock_pause = true;
                    std::cout << "unpaused and locked pause" << std::endl;
                }
                break;
            }
            break;
        case SDL_KEYUP:
        {
            switch (event.key.keysym.sym)
            {
            case SDLK_SPACE:
            {
                if (lock_pause)
                {
                    lock_pause = false;
                    std::cout << "lock released" << std::endl;
                }
                break;
            }
            case SDLK_1:
            {
                //homo
                windowSurface = SDL_GetWindowSurface(window);
                SDL_FillRect(windowSurface, NULL, 0x000000);
                int w;
                int h;
                SDL_GetWindowSize(window, &w, &h);
                int rows;
                int columns;
                std::cout << "how many rows?";
                std::cin >> rows;
                std::cout << "how many columns?";
                std::cin >> columns;

                double w_gap = w / (columns + 1);
                double h_gap = h / (rows + 1);

                for (double y = h_gap; y <= h_gap * rows; y += h_gap)
                    for (double x = w_gap; x <= w_gap * columns; x += w_gap)
                        setPixelRandomColor(renderer, windowSurface, floor((int)x), floor((int)y));

                auto texture = SDL_CreateTextureFromSurface(renderer, windowSurface);
                SDL_RenderCopy(renderer, texture, NULL, NULL);
                SDL_RenderPresent(renderer);
                SDL_DestroyTexture(texture);
                break;
            }
            case SDLK_2:
            {
                //range
                windowSurface = SDL_GetWindowSurface(window);
                SDL_FillRect(windowSurface, NULL, 0x000000);
                int w;
                int h;
                SDL_GetWindowSize(window, &w, &h);
                int seed_count;
                std::cout << "how many seeds?";
                std::cin >> seed_count;


                for (int i = seed_count; i > 0; i--)
                {
                    setPixelRandomColor(renderer, windowSurface, random_coord(w), random_coord(h));
                }

                auto texture = SDL_CreateTextureFromSurface(renderer, windowSurface);
                SDL_RenderCopy(renderer, texture, NULL, NULL);
                SDL_RenderPresent(renderer);
                SDL_DestroyTexture(texture);
                break;
            }
            case SDLK_3:
            {
                //rand
                windowSurface = SDL_GetWindowSurface(window);
                std::vector<Point> points;
                SDL_FillRect(windowSurface, NULL, 0x000000);
                int w;
                int h;
                SDL_GetWindowSize(window, &w, &h);
                int seed_count;
                int range;
                std::cout << "how many seeds?";
                std::cin >> seed_count;
                std::cout << "range?";
                std::cin >> range;

                int repeated_iteration = 0;
                do {
                    int x = random_coord(w);
                    int y = random_coord(h);
                    Point p(x, y);
                    if (none_of(points.begin(), points.end(), [p, range](Point pp) { return pp.isInRadius(&p, range); }))
                    {
                        setPixelRandomColor(renderer, windowSurface, x, y);
                        repeated_iteration = 0;
                        seed_count--;
                        points.push_back(p);
                    }
                    else repeated_iteration++;
                } while (seed_count > 0 && repeated_iteration < 500);
                if (seed_count > 0) std::cout << "failed to draw " << seed_count << " seeds" << std::endl;
                auto texture = SDL_CreateTextureFromSurface(renderer, windowSurface);
                SDL_RenderCopy(renderer, texture, NULL, NULL);
                SDL_RenderPresent(renderer);
                SDL_DestroyTexture(texture);
                break;
            }
            case SDLK_4:
            {
                int new_window_w;
                int new_window_h;
                std::cout << "Width: ";
                std::cin >> new_window_w;
                std::cout << "Height: ";
                std::cin >> new_window_h;
                screen_w = new_window_w;
                screen_h = new_window_h;
                SDL_SetWindowSize(window, new_window_w, new_window_h);
                break;
            }
            case SDLK_5:
            {
                int input;
                std::cout << "0 - cancel 1 - no period 2 - period\n";
                std::cin >> input;
                if (input == 0) break;
                else if (input == 1) periodicalness = false;
                else if (input == 2) periodicalness = true;
                break;
            }
            case SDLK_6:
            {
                int input;
                std::cout << "0: vonNeumann\n1:Moore\n2: pentagonal left\n3: pentagonal right\n4: pentagonal up\n5: pentagonal down\n6: pentagonal random\n7:hexagonal left\n8: hexagonal right\n9:hexagonal random\n10: cancel";
                std::cin >> input;
                if (input == 0) break;
                switch (input)
                {
                case 0:
                    nucleation_mode = VONNEUMAN;
                    break;
                case 1:
                    nucleation_mode = MOORE;
                    break;
                case 2:
                    nucleation_mode = PENTALEFT;
                    break;
                case 3:
                    nucleation_mode = PENTARIGHT;
                    break;
                case 4:
                    nucleation_mode = PENTAUP;
                    break;
                case 5:
                    nucleation_mode = PENTADOWN;
                    break;
                case 6:
                    nucleation_mode = PENTARANDOM;
                    break;
                case 7:
                    nucleation_mode = HEXALEFT;
                    break;
                case 8:
                    nucleation_mode = HEXARIGHT;
                    break;
                case 9:
                    nucleation_mode = HEXARANDOM;
                    break;
                default:
                    break;
                }
                break;
            }
            case SDLK_7:
            {
                std::cout << "7" << std::endl;
                if (!montecarlo)
                {
                    char choise;
                    std::cout << "do you want to turn on monte carlo? y/n " << std::endl;
                    std::cin >> choise;
                    if (choise != 'y') break;
                    else
                    {
                        std::cout << "montecarlo = 1\n";
                        borders.clear();
                        //searchForBorders(windowSurface);
                        montecarlo = true;
                    }
                }
                else
                {

                    char choise;
                    std::cout << "do you want to turn off monte carlo? y/n " << std::endl;
                    std::cin >> choise;
                    if (choise == 'y')
                    {
                        montecarlo = false;

                        std::cout << "montecarlo = 0\n";
                    }
                }
            }

            case SDLK_8:
            {
                std::cout << "7" << std::endl;
                if (!montecarlo)
                {
                    char choise;
                    std::cout << "do you want to turn on recrystalisation? y/n " << std::endl;
                    std::cin >> choise;
                    if (choise != 'y') break;
                    else
                    {
                        std::cout << "recrystalisation = 1\n";
                        //searchForBorders(windowSurface);
                        recrystalisation = true;
                    }
                }
            }
            }
            break;
        }
        case SDL_MOUSEBUTTONUP:
        {
            if (click != true)
                click = true;
            else break;
            //if (pause) break;
            int x;
            int y;
            bool loop;
            SDL_GetMouseState(&x, &y);
            setPixelRandomColor(renderer, windowSurface, x, y);
            auto texture = SDL_CreateTextureFromSurface(renderer, windowSurface);
            SDL_RenderCopy(renderer, texture, NULL, NULL);
            SDL_RenderPresent(renderer);
            SDL_DestroyTexture(texture);
            break;
        }
        case SDL_QUIT:
            quit = true;
            break;
        default:
        {
            if (pause)
                break;
            if (recrystalisation)
            {
                if(recrystalise(window, renderer)) 
                    recrystalisation = false;
                break;
            }
            if (montecarlo)
            {
                monteCarlo(windowSurface, 0);
                auto texture = SDL_CreateTextureFromSurface(renderer, windowSurface);
                SDL_RenderCopy(renderer, texture, NULL, NULL);
                SDL_RenderPresent(renderer);
                SDL_DestroyTexture(texture);
                break;
            }
            int w;
            int h;
            SDL_GetWindowSize(window, &w, &h);
            SDL_Surface* oldSurface = SDL_CreateRGBSurface(0, w, h, 32, 0, 0, 0, 0);
            SDL_BlitSurface(windowSurface, NULL, oldSurface, NULL);

            for (int i = 0; i < w; i++)
                for (int j = 0; j < h; j++)
                {
                    int pxColor = probePixel(oldSurface, i, j, false);
                    if (pxColor != 0) continue;
                    std::vector<std::pair<int, int>> neighbors;
                    switch (nucleation_mode)
                    {
                    case VONNEUMAN:

                        break;
                    }

                    switch (nucleation_mode)
                    {
                    case VONNEUMAN:
                        vonNeumann(oldSurface, i, j, &neighbors, periodicalness);
                        break;
                    case MOORE:
                        moore(oldSurface, i, j, &neighbors, periodicalness);
                        break;
                    case PENTALEFT:
                        pentaLeft(oldSurface, i, j, &neighbors, periodicalness);
                        break;
                    case PENTARIGHT:
                        pentaRight(oldSurface, i, j, &neighbors, periodicalness);
                        break;
                    case PENTAUP:
                        pentaUp(oldSurface, i, j, &neighbors, periodicalness);
                        break;
                    case PENTADOWN:
                        pentaDown(oldSurface, i, j, &neighbors, periodicalness);
                        break;
                    case PENTARANDOM:
                    {
                        int n = random_coord(4);
                        switch (n)
                        {
                        case 0:
                            pentaLeft(oldSurface, i, j, &neighbors, periodicalness);
                            break;
                        case 1:
                            pentaRight(oldSurface, i, j, &neighbors, periodicalness);
                            break;
                        case 2:
                            pentaUp(oldSurface, i, j, &neighbors, periodicalness);
                            break;
                        case 3:
                            pentaDown(oldSurface, i, j, &neighbors, periodicalness);
                            break;
                        }
                        break;
                    }
                    case HEXALEFT:
                        hexaLeft(oldSurface, i, j, &neighbors, periodicalness);
                        break;
                    case HEXARIGHT:
                        hexaRight(oldSurface, i, j, &neighbors, periodicalness);
                        break;
                    case HEXARANDOM:
                    {
                        int n = random_coord(2);
                        switch (n)
                        {
                        case 0:
                            hexaLeft(oldSurface, i, j, &neighbors, periodicalness);
                            break;
                        case 1:
                            hexaRight(oldSurface, i, j, &neighbors, periodicalness);
                            break;
                        }
                        break;
                    }
                    default:
                        break;
                    }
                    //if (montecarlo) break;
                    std::vector<int> colorToSet;
                    for (std::pair<int, int> n : neighbors)
                        if (std::all_of(neighbors.begin(), neighbors.end(), [colorToSet, n](std::pair<int, int> nn) { if (n.first == nn.first || n.second > nn.second) return true; else return false; }))
                        {
                            setPixel(windowSurface, i, j, n.first);
                            break;
                        }
                        else
                        {
                            int max = 0;
                            std::for_each(neighbors.begin(), neighbors.end(),
                                [&colorToSet, max](std::pair<int, int> p)
                                {
                                    if (p.second > max)
                                    {
                                        colorToSet.clear();
                                        colorToSet.push_back(p.first);
                                    }
                                    else if (p.second == max)
                                        colorToSet.push_back(p.first);
                                }
                            );
                        }
                    if (colorToSet.size() > 0) setPixel(windowSurface, i, j, colorToSet[random_coord(colorToSet.size())]);
                }
            auto texture = SDL_CreateTextureFromSurface(renderer, windowSurface);

            SDL_RenderCopy(renderer, texture, NULL, NULL);
            SDL_RenderPresent(renderer);
            SDL_DestroyTexture(texture);
            SDL_FreeSurface(oldSurface);
            break;
        }
        }
    }
    //end
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

}
int main(int argc, char** argv)
{
    mainLoop();

    return 0;
}