/*
 * dcdio.h
 *
 *  Created on: Feb 27, 2009
 *      Author: zhmurov
 */

#pragma once

#include <stdio.h>

struct DCD{
    FILE* file;
    int N, NFILE, NPRIV, NSAVC;
    double DELTA;
    float *X, *Y, *Z;

    void open_read(const char* filename);
    void open_write(const char* filename);
    void open_append(const char* filename);
    void close();

    void write_header() const;
    void read_header();
    void write_frame(float *X, float *Y, float *Z) const;
    void write_frame() const { write_frame(this->X, this->Y, this->Z); }
    int read_frame(float *X, float *Y, float *Z);
    int read_frame() { return read_frame(this->X, this->Y, this->Z); }

    void goto_frame(int frame);
    void goto_header();

    void allocate();
    void deallocate();

    DCD() : X(NULL), Y(NULL), Z(NULL) { }
    ~DCD() { deallocate(); }
};

