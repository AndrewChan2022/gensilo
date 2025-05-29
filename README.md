# gensilo

## dependency

download silo+hdf5 from 

https://visit-dav.github.io/visit-website/releases-as-tables/

the dev link:

https://github.com/visit-dav/visit/releases/download/v3.4.1/visit_windowsdev_3.4.1.zip

unzip and get sili+hdf5


## build

<!-- make sure Silo/bin before visit location -->

change you lib location:

```bash

if (MSVC)
    set(CMAKE_PREFIX_PATH  
        "E:/local/silo4.10.3" 
        "E:/local/hdf5/1.8.19" 
        ${CMAKE_PREFIX_PATH} 
    )
endif()


```


```bash

mkdir build
cd build
cmake ..


```

