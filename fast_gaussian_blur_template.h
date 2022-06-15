// Copyright (C) 2017-2022 Basile Fraboni
// Copyright (C) 2014 Ivan Kutskir (for the original fast blur implmentation)
// All Rights Reserved
// You may use, distribute and modify this code under the
// terms of the MIT license. For further details please refer 
// to : https://mit-license.org/
//
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
typedef unsigned char uchar;
//!
//! \file fast_gaussian_blur_template.h
//! \author Basile Fraboni
//! \date 2017 - 2022
//!
//! \brief This contains a C++ implementation of a fast Gaussian blur algorithm in linear time.
//!
//! The image buffer is supposed to be of size w * h * c, h its height, with w its width, 
//! and c its number of channels.
//! The default implementation only supports up to 4 channels images, but one can easily add support for any number of channels
//! using either specific template cases or a generic function that takes the number of channels as an explicit parameter.
//! This implementation is focused on learning and readability more than on performance.
//! The fast blur algorithm is performed with several box blur passes over an image.
//! The filter converges towards a true Gaussian blur after several passes (thanks TCL). In practice,
//! three passes are sufficient for good quality results.
//! For further details please refer to:
//!     - http://blog.ivank.net/fastest-gaussian-blur.html
//!     - https://www.peterkovesi.com/papers/FastGaussianSmoothing.pdf
//!     - https://github.com/bfraboni/FastGaussianBlur
//!
//! **Note:** The fast gaussian blur algorithm is not accurate on image boundaries. 
//! It performs a diffusion of the signal with several independant passes, each pass depending 
//! of the preceding one. Some of the diffused signal is lost near borders and results in a slight 
//! loss of accuracy for next pass. This problem can be solved by increasing the image support of 
//! half the box kernel extent at each pass of the algorithm. The added padding would in this case 
//! capture the diffusion and make the next pass accurate. 
//! On contrary true Gaussian blur does not suffer this problem since the whole diffusion process 
//! is performed in one pass only.
//! The extra padding is not performed in this implementation, however we provide and discuss several border
//! policies resulting in dfferent approximations and accuracies.  
//! 


//!
//! \fn void std_to_box(float boxes[], float sigma, int n)  
//!
//! \brief this function converts the standard deviation of 
//! Gaussian blur into dimensions of boxes for box blur. For 
//! further details please refer to :
//! https://www.peterkovesi.com/matlabfns/#integral
//! https://www.peterkovesi.com/papers/FastGaussianSmoothing.pdf
//!
//! \param[out] boxes   boxes dimensions
//! \param[in] sigma    Gaussian standard deviation
//! \param[in] n        number of boxes
//!
void std_to_box(int boxes[], float sigma, int n)  
{
    // ideal filter width
    float wi = sqrt((12*sigma*sigma/n)+1); 
    int wl = wi; // no need std::floor  
    if(wl%2==0) wl--;
    int wu = wl+2;
                
    float mi = (12*sigma*sigma - n*wl*wl - 4*n*wl - 3*n)/(-4*wl - 4);
    int m = mi+0.5f; // avoid std::round by adding 0.5f and cast to integer type

    for(int i=0; i<n; i++) 
        boxes[i] = ((i < m ? wl : wu) - 1) / 2;
}

//!
//! \fn void horizontal_blur(int * in, int * out, int w,int h, int r)   
//!
//! \brief this function performs the horizontal blur pass for box blur. 
//!
//! \param[in,out] in       source channel
//! \param[in,out] out      target channel
//! \param[in] w            image width
//! \param[in] h            image height
//! \param[in] r            box dimension
//!
void horizontal_blur(uchar * in, uchar * out, int w, int h, int r)
{
    float iarr = 1.f / (r+r+1);
    for(int i=0; i<h; i++) 
    {
        int ti = i*w; 
        int li = ti;
        int ri = ti+r;
        int fv = in[ti];
        int lv = in[ti+w-1];
        int val = (r+1)*fv;

        for(int j=0;   j<r;   j++) { 
            val += in[ti+j]; 
        }
        
        for(int j=0  ; j<=r ; j++) { 
            val += in[ri++] - fv; 
            out[ti++] = round(val*iarr); 
        }
        
        for(int j=r+1; j<w-r; j++) { 
            val += in[ri++] - in[li++]; 
            out[ti++] = round(val*iarr); 
        }
        
        for(int j=w-r; j<w  ; j++) { 
            val += lv - in[li++]; 
            out[ti++] = round(val*iarr); 
        }
    }
}

//!
//! \fn void total_blur(int * in, int * out, int w, int h, int r)   
//!
//! \brief this function performs the total blur pass for box blur. 
//!
//! \param[in,out] in       source channel
//! \param[in,out] out      target channel
//! \param[in] w            image width
//! \param[in] h            image height
//! \param[in] r            box dimension
//!
void total_blur(uchar * in, uchar * out, int w, int h, int r) 
{
    float iarr = 1.f / (r+r+1);
    for(int i=0; i<w; i++) 
    {
        int ti = i;
        int li = ti;
        int ri = ti+r*w;
        int fv = in[ti];
        int lv = in[ti+w*(h-1)];
        int val = (r+1)*fv;

        for(int j=0;   j<r;   j++) { 
            val += in[ti+j*w]; 
        }
        
        for(int j=0  ; j<=r ; j++) { 
            val += in[ri] - fv; 
            out[ti] = round(val*iarr);
            ri+=w; ti+=w;
        }
        
        for(int j=r+1; j<h-r; j++) { 
            val += in[ri] - in[li]; 
            out[ti] = round(val*iarr); 
            li+=w; 
            ri+=w; 
            ti+=w; 
        }
        
        for(int j=h-r; j<h  ; j++) { 
            val += lv - in[li]; 
            out[ti] = round(val*iarr); 
            li+=w; 
            ti+=w; 
        }
    }
}

void swap(void  *v1, void  *v2, size_t size)  
{
#if 0
    void *temp = malloc(size);
    assert(temp != NULL);
    memcpy(temp, v1, size);
    memcpy(v1, v2, size);
    memcpy(v2, temp, size);
    free(temp);
#endif
    uchar *p1 = (uchar *)v1;
    uchar *p2 = (uchar *)v2;
    size_t i;
    for (i = 0; i < size; i++) {
        uchar tmp = p1[i];
        p1[i] = p2[i];
        p2[i] = tmp;
    }
}

void box_blur(uchar * in, uchar * out, size_t size, int w, int h, int r) 
{
    swap(in, out, size);
    horizontal_blur(out, in, w, h, r);
    total_blur(in, out, w, h, r);
    // Note to myself : 
    // here we could go anisotropic with different radiis rx,ry in HBlur and TBlur
}

void fast_gaussian_blur_int(uchar * in, uchar * out, size_t size, int w, int h, float sigma) 
{
    // sigma conversion to box dimensions
    int boxes[3];
    std_to_box(boxes, sigma, 3);
    box_blur(in, out, size, w, h, boxes[0]);
    box_blur(out, in, size, w, h, boxes[1]);
    box_blur(in, out, size, w, h, boxes[2]);
}

//!
//! \fn void horizontal_blur_rgb(int * in, int * out, int w, int h, int c, int r)   
//!
//! \brief this function performs the horizontal blur pass for box blur. 
//!
//! \param[in,out] in       source channel
//! \param[in,out] out      target channel
//! \param[in] w            image width
//! \param[in] h            image height
//! \param[in] c            image channels
//! \param[in] r            box dimension
//!
void horizontal_blur_rgb(uchar * in, uchar * out, int w, int h, int c, int r) 
{
    // radius range on either side of a pixel + the pixel itself
    float iarr = 1.f / (r+r+1);
    #pragma omp parallel for
    for(int i=0; i<h; i++) 
    {
        int ti = i*w; 
        int li = ti;  
        int ri = ti+r;
        //1, 0, 61765110
        int tmp_fv[3] = {1,0,61765110};
        int fv[3] = { in[ti*c+0], in[ti*c+1], in[ti*c+2] };                  
        int lv[3] = { in[(ti+w-1)*c+0], in[(ti+w-1)*c+1], in[(ti+w-1)*c+2] };
        int val[3] = { (r+1)*fv[0], (r+1)*fv[1], (r+1)*fv[2] };              

        for(int j=0; j<r; j++) 
        { 
            val[0] += in[(ti+j)*c+0]; 
            val[1] += in[(ti+j)*c+1]; 
            val[2] += in[(ti+j)*c+2]; 
        }

        for(int j=0; j<=r; j++, ri++, ti++) 
        { 
            val[0] += in[ri*c+0] - fv[0]; 
            val[1] += in[ri*c+1] - fv[1]; 
            val[2] += in[ri*c+2] - fv[2]; 
            out[ti*c+0] = round(val[0]*iarr); 
            out[ti*c+1] = round(val[1]*iarr); 
            out[ti*c+2] = round(val[2]*iarr); 
        }

        for(int j=r+1; j<w-r; j++, ri++, ti++, li++) 
        { 
            val[0] += in[ri*c+0] - in[li*c+0]; 
            val[1] += in[ri*c+1] - in[li*c+1]; 
            val[2] += in[ri*c+2] - in[li*c+2]; 
            out[ti*c+0] = round(val[0]*iarr); 
            out[ti*c+1] = round(val[1]*iarr); 
            out[ti*c+2] = round(val[2]*iarr); 
        }

        for(int j=w-r; j<w; j++, ti++, li++) 
        { 
            val[0] += lv[0] - in[li*c+0]; 
            val[1] += lv[1] - in[li*c+1]; 
            val[2] += lv[2] - in[li*c+2]; 
            out[ti*c+0] = round(val[0]*iarr); 
            out[ti*c+1] = round(val[1]*iarr); 
            out[ti*c+2] = round(val[2]*iarr); 
        }
    }
}

//!
//! \fn void total_blur_rgb(int * in, int * out, int w, int h, int c, int r)   
//!
//! \brief this function performs the total blur pass for box blur. 
//!
//! \param[in,out] in       source channel
//! \param[in,out] out      target channel
//! \param[in] w            image width
//! \param[in] h            image height
//! \param[in] c            image channels
//! \param[in] r            box dimension
//!
void total_blur_rgb(uchar * in, uchar * out, int w, int h, int c, int r) 
{
    // radius range on either side of a pixel + the pixel itself
    float iarr = 1.f / (r+r+1);
    for(int i=0; i<w; i++) 
    {
        int ti = i;
        int li = ti;
        int ri = ti+r*w;

        int fv[3] = {in[ti*c+0], in[ti*c+1], in[ti*c+2] };
        int lv[3] = {in[(ti+w*(h-1))*c+0], in[(ti+w*(h-1))*c+1], in[(ti+w*(h-1))*c+2] };
        int val[3] = {(r+1)*fv[0], (r+1)*fv[1], (r+1)*fv[2] };

        for(int j=0; j<r; j++) 
        { 
            val[0] += in[(ti+j*w)*c+0]; 
            val[1] += in[(ti+j*w)*c+1]; 
            val[2] += in[(ti+j*w)*c+2]; 
        }

        for(int j=0; j<=r; j++, ri+=w, ti+=w) 
        { 
            val[0] += in[ri*c+0] - fv[0]; 
            val[1] += in[ri*c+1] - fv[1]; 
            val[2] += in[ri*c+2] - fv[2]; 
            out[ti*c+0] = round(val[0]*iarr); 
            out[ti*c+1] = round(val[1]*iarr); 
            out[ti*c+2] = round(val[2]*iarr); 
        }

        for(int j=r+1; j<h-r; j++, ri+=w, ti+=w, li+=w) 
        { 
            val[0] += in[ri*c+0] - in[li*c+0]; 
            val[1] += in[ri*c+1] - in[li*c+1]; 
            val[2] += in[ri*c+2] - in[li*c+2]; 
            out[ti*c+0] = round(val[0]*iarr); 
            out[ti*c+1] = round(val[1]*iarr); 
            out[ti*c+2] = round(val[2]*iarr); 
        }
        
        for(int j=h-r; j<h; j++, ti+=w, li+=w) 
        { 
            val[0] += lv[0] - in[li*c+0]; 
            val[1] += lv[1] - in[li*c+1]; 
            val[2] += lv[2] - in[li*c+2]; 
            out[ti*c+0] = round(val[0]*iarr); 
            out[ti*c+1] = round(val[1]*iarr); 
            out[ti*c+2] = round(val[2]*iarr); 
        }
    }
}

//!
//! \fn void box_blur_rgb(int * in, int * out, int w, int h, int c, int r)   
//!
//! \brief this function performs a box blur pass. 
//!
//! \param[in,out] in       source channel
//! \param[in,out] out      target channel
//! \param[in] w            image width
//! \param[in] h            image height
//! \param[in] c            image channels
//! \param[in] r            box dimension
//!
void box_blur_rgb(uchar * in, uchar * out, size_t size, int w, int h, int c, int r) 
{
    swap(in, out, size);
    horizontal_blur_rgb(out, in, w, h, c, r);
    total_blur_rgb(in, out, w, h, c, r);
    // Note to myself : 
    // here we could go anisotropic with different radiis rx,ry in HBlur and TBlur
}

//!
//! \fn void fast_gaussian_blur_rgb(int * in, int * out, int w, int h, int c, float sigma)   
//!
//! \brief this function performs a fast Gaussian blur. Applying several
//! times box blur tends towards a true Gaussian blur. Three passes are sufficient
//! for good results. For further details please refer to :  
//! http://blog.ivank.net/fastest-gaussian-blur.html
//!
//! \param[in,out] in       source channel
//! \param[in,out] out      target channel
//! \param[in] w            image width
//! \param[in] h            image height
//! \param[in] c            image channels
//! \param[in] sigma        gaussian std dev
//!
void fast_gaussian_blur_rgb(uchar * in, uchar * out, size_t size, int w, int h, int c, float sigma) 
{
    // sigma conversion to box dimensions
    int boxes[3];
    std_to_box(boxes, sigma, 3);
    box_blur_rgb(in, out, size, w, h, c, boxes[0]);
    box_blur_rgb(out, in, size, w, h, c, boxes[1]);
    box_blur_rgb(in, out, size, w, h, c, boxes[2]);
}

void horizontal_blur_t (uchar * in, uchar * out, int w, int h, int c, int r) 
{
    float iarr = 1.f / (r+r+1);
    #pragma omp parallel for
    for(int i=0; i<h; i++) 
    {
        int ti = i*w, li = ti, ri = ti+r;
        //float fv[c] , lv[c], val[c];
        float  val[c];

        for(int ch = 0; ch < c; ++ch)
        {            
            fv[ch]  =  0;
            lv[ch]  =  0;
            val[ch] =  0;
        }

        // initial acucmulation
        for(int j=0; j<r; j++) 
        for(int ch = 0; ch < c; ++ch)
        {
            val[ch] += in[(ti+j)*c+ch]; 
        }

        // left border - filter kernel is incomplete
        for(int j=0; j<=r; j++, ri++, ti++) 
        for(int ch = 0; ch < c; ++ch)
        { 
            val[ch] +=  in[ri*c+ch]; 
            out[ti*c+ch] = val[ch]/(r+j+1); 
        }

        // center of the image - filter kernel is complete
        for(int j=r+1; j<w-r; j++, ri++, ti++, li++) 
        for(int ch = 0; ch < c; ++ch)
        { 
            val[ch] += in[ri*c+ch] - in[li*c+ch]; 
            out[ti*c+ch] = val[ch]*iarr;
        }

        // right border - filter kernel is incomplete
        for(int j=w-r; j<w; j++, ti++, li++) 
        for(int ch = 0; ch < c; ++ch)
        { 
            val[ch] += (-in[li*c+ch]); 
            out[ti*c+ch] =  (val[ch]/(r+w-j)); 
        }
    }
}

#define MIN(a,b) (((a) < (b)) ? (a) : (b))  
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

void flip_block(uchar * in, uchar * out, int w,  int h, int c)
{
    const int block = 256/c;
    #pragma omp parallel for collapse(2)
    for(int x= 0; x < w; x+= block)
    for(int y= 0; y < h; y+= block)
    {
        const uchar * p = in + y*w*c + x*c;
        uchar * q = out + y*c + x*h*c;
        int tmp1 = (x + block);
        const int blockx= MIN(w, tmp1) - x;
        int tmp2 = (y + block);
        const int blocky= MIN(h, tmp2) - y;
        for(int xx= 0; xx < blockx; xx++)
        {
            for(int yy= 0; yy < blocky; yy++)
            {
                for(int k= 0; k < c; k++)
                    q[k]= p[k];
                p+= w*c;
                q+= c;                    
            }
            p+= -blocky*w*c + c;
            q+= -blocky*c + h*c;
        }
    }
}

#if 0
void fast_gaussian_blur_t(uchar * in, uchar * out, size_t size, int w, int h, int c, float sigma) 
{
    // compute box kernel sizes
    int boxes[3];

    struct timeval tv_start_1, tv_end_1;
    struct timeval tv_start_2, tv_end_2;
    struct timeval tv_start_3, tv_end_3;
    struct timeval tv_start_4, tv_end_4;

    // stats
    gettimeofday(&tv_start_1, NULL);
    std_to_box(boxes, sigma, 3);
    // stats
    gettimeofday(&tv_end_1, NULL);

    float elapsed_1 = ((tv_end_1.tv_sec - tv_start_1.tv_sec) +  1e-6 * (tv_end_1.tv_usec - tv_start_1.tv_usec)) * 1000.00;
    printf("c time 1 %.4fms\n", elapsed_1);


        // stats
    gettimeofday(&tv_start_2, NULL);
    // perform 3 horizontal blur passes
    horizontal_blur_t(in, out, w, h, c, boxes[0]);
    horizontal_blur_t(out, in, w, h, c, boxes[1]);
    horizontal_blur_t(in, out, w, h, c, boxes[2]);
    // stats
    gettimeofday(&tv_end_2, NULL);

    float elapsed_2 = ((tv_end_2.tv_sec - tv_start_2.tv_sec) +  1e-6 * (tv_end_2.tv_usec - tv_start_2.tv_usec)) * 1000.00;
    printf("c time 2 %.4fms\n", elapsed_2);

            // stats
    gettimeofday(&tv_start_3, NULL);
    // flip buffer
    flip_block(out, in, w, h, c);
        // stats
    gettimeofday(&tv_end_3, NULL);
    
    float elapsed_3 = ((tv_end_3.tv_sec - tv_start_3.tv_sec) +  1e-6 * (tv_end_3.tv_usec - tv_start_3.tv_usec)) * 1000.00;
    printf("c time 3 %.4fms\n", elapsed_3);
    
    // perform 3 horizontal blur passes on flipped image
    horizontal_blur_t(in, out, h, w, c, boxes[0]);
    horizontal_blur_t(out, in, h, w, c, boxes[1]);
    horizontal_blur_t(in, out, h, w, c, boxes[2]);
    
    // flip buffer
    flip_block(out, in, h, w, c);

    // stats
    gettimeofday(&tv_start_4, NULL);
    // swap pointers to get result in the ouput buffer 
    swap(in, out, size);
    gettimeofday(&tv_end_4, NULL);

    float elapsed_4 = ((tv_end_4.tv_sec - tv_start_4.tv_sec) +  1e-6 * (tv_end_4.tv_usec - tv_start_4.tv_usec)) * 1000.00;
    printf("c time 4 %.4fms\n", elapsed_4);

}
#endif

void fast_gaussian_blur_t(uchar * in, uchar * out, size_t size, int w, int h, int c, float sigma) 
{
    // compute box kernel sizes
    int boxes[3];

    struct timeval tv_start_1, tv_end_1;
    struct timeval tv_start_2, tv_end_2;
    struct timeval tv_start_3, tv_end_3;
    struct timeval tv_start_4, tv_end_4;

    std_to_box(boxes, sigma, 3);

    // perform 3 horizontal blur passes
    horizontal_blur_t(in, out, w, h, c, boxes[0]);
    horizontal_blur_t(out, in, w, h, c, boxes[1]);
    horizontal_blur_t(in, out, w, h, c, boxes[2]);

    // flip buffer
    flip_block(out, in, w, h, c);

    // perform 3 horizontal blur passes on flipped image
    horizontal_blur_t(in, out, h, w, c, boxes[0]);
    horizontal_blur_t(out, in, h, w, c, boxes[1]);
    horizontal_blur_t(in, out, h, w, c, boxes[2]);
    
    // flip buffer
    flip_block(out, in, h, w, c);

    // swap pointers to get result in the ouput buffer 
    swap(in, out, size);
}
