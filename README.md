# Autoimmune Landscape project


## Data Preparation Instructions

Before running the code in this repository, follow the steps below to set up the necessary datasets.

### 1. Download the Dataset Folder

Download the `downloaded_datasets` folder from the following Google Drive link:

https://drive.google.com/drive/folders/1ykQ6lj9FOASaWdjL3ia6ePKk4HaHsODD?usp=sharing

Then place it in the following location within your project directory:

```
AutoimmuneLandscape/Data_Collection/
```

Your folder structure should look like this:

```
AutoimmuneLandscape/
└── Data_Collection/
    └── downloaded_datasets/
        ├── crohn/
        │   ├── downloaded_microarray/
        │   └── downloaded_RNAseq/
        └── ...
```

Only `.xz` compressed files are included to minimize storage usage.

### 2. Uncompress All `.xz` Files

To uncompress all `.xz` files while keeping the original compressed files, run the following command from the root of the project (`AutoimmuneLandscape/`):

```bash
find Data_Collection/downloaded_datasets/ -type f -name "*.xz" -exec unxz -k {} \;
```

This will decompress all `.xz` files across the directory tree.

