#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include<sys/time.h>

// image io
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

// fast blur
#include "fast_gaussian_blur_template.h"

typedef unsigned char uchar;

// #define MIN(a,b) (((a) < (b)) ? (a) : (b))  
// #define MAX(a,b) (((a) > (b)) ? (a) : (b))


// #define USE_FLOAT
typedef int16_t lv_coord_t;
typedef struct {
    lv_coord_t x1;
    lv_coord_t y1;
    lv_coord_t x2;
    lv_coord_t y2;
} lv_area_t;

int main(int argc, char * argv[])
{   
    // helper
    if( argc < 4 )
    {
        printf("%s [input] [output] [sigma] [passes - optional]\n", argv[0]);
        exit(1);
    }
    
    // load image
    int width, height, channels;
    uchar * image_data = stbi_load(argv[1], &width, &height, &channels, 0);
    printf("Source image: %s %dx%d (%d)\n", argv[1], width, height, channels);

    // read parameters
    float sigma = atof(argv[3]);
    int passes = argc > 4 ? atoi(argv[4]) : 3;
    
    // temporary data
    size_t size = width * height * channels;
    size_t i;


    //int new_image[size];
    //int old_image[size];

    uchar * new_image = (uchar *)malloc(size * sizeof(uchar));
    uchar * old_image = (uchar *)malloc(size * sizeof(uchar));
    
    memset(new_image, 0, size*sizeof(uchar));
    memset(old_image, 0, size*sizeof(uchar));

    FILE  *fd = fopen("blur_begin","w+");
    // channels copy r,g,b
    for(i = 0; i < size; i++) {
        old_image[i] = image_data[i];
        fprintf(fd, "%02x ", old_image[i]); 
       if ( ( i % 1920 == 0) && (i != 0) ) {
           fprintf(fd, "\n");
       }
    }
    fclose(fd);

    struct timeval tv_start, tv_end;
    
    // stats
    gettimeofday(&tv_start, NULL);

    // perform gaussian blur
    // note: the implementation can work on any buffer types (uint8, uint16, uint32, int, float, double)
    // note: both old and new buffer are modified
    //fast_gaussian_blur(old_image, new_image, width, height, channels, sigma, passes);
    //fast_gaussian_blur_int(old_image, new_image, size, width, height, sigma);
    //fast_gaussian_blur_rgb(old_image, new_image, size, width, height, channels, sigma);


    int image_size =  480 * 480 * channels * sizeof(float);

    uint8_t  *old_image_buf = malloc(image_size);
    if(old_image_buf == NULL) return;

    uint8_t  *new_image_buf = malloc(image_size);
    if(new_image_buf == NULL) return;

    lv_area_t a = {24, 311, 24 + 208, 0 + 311+145};

    a.x1 =  a.x1 * channels; 
    a.x2 = a.x1 + a.x2  * channels;

    int loop = 0; 
    for (int h = a.y1; h < a.y2; h++) {
        for (int w = a.x1; w < a.x2; w++) {
            int index = h * 480 * channels + w; /*Calculate the position of the first pixel value of the object to be blurred*/
            old_image_buf[loop] = old_image[index];
            loop++;
        }
    }
    
    fast_gaussian_blur_t(old_image_buf, new_image_buf, size, 208, 145, channels, sigma);
    
    // stats
    gettimeofday(&tv_end, NULL);

    float elapsed = ((tv_end.tv_sec - tv_start.tv_sec) +  1e-6 * (tv_end.tv_usec - tv_start.tv_usec)) * 1000.00;
    printf("Time %.4fms\n", elapsed);

    loop = 0;
    for (int h = a.y1; h < a.y2; h++) {
        for (int w = a.x1; w < a.x2; w++) {
            int index = h * 480 * channels + w;
            old_image[index] = new_image_buf[loop];
            loop++;
        }
    }



    FILE  *fd_1 = fopen("blur_end","w+");

    // convert result
    for(i = 0; i < size; i++) {
        image_data[i] = (uchar)old_image[i];
        
        fprintf(fd_1, "%02x ", old_image[i]);
       
        if ( ( i % 1920 == 0) && (i != 0) ) {
            fprintf(fd_1, "\n");
        }
    }
    fclose(fd_1);



    // save image
    char image_type[256];
    memset(image_type, 0, sizeof(image_type));
    size_t image_len = strlen(argv[2]);
    strncpy(image_type, argv[2] + image_len - 3, 3);
    
    if (strcmp(image_type, "bmp") == 0) {
        stbi_write_bmp(argv[2], width, height, channels, image_data);
    } else if (strcmp(image_type, "jpg") == 0) {
        stbi_write_jpg(argv[2], width, height, channels, image_data, 90);
    } else {

        if (strcmp(image_type, "png") != 0) {
            printf("Image format '%s' not supported, writing default png\n", image_type);
            stbi_write_png("default.png", width, height, channels, image_data, channels*width);
        } else {
            stbi_write_png(argv[2], width, height, channels, image_data, channels*width);
        }
    }  

    free(new_image);
    free(old_image);
    
    return 0;
}
