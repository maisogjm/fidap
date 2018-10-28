struct analyze_struct
{
	struct header_key
	{
		int sizeof_hdr;
		char data_type[10];
		char db_name[18];
		int extents;
		short int session_error;
		char regular;
		char hkey_un0;
	} hk;
	struct image_dimension                  /*      image_dimension  */
	{                                   /* off + size*/
	short int dim[8];               /* 0 + 16    */
        char vox_units[4];              /* 16 + 4    */
        char cal_units[8];              /* 20 + 4    */
	short int unused1;              /* 24 + 2    */
        short int datatype;             /* 30 + 2    */
	short int bitpix;               /* 32 + 2    */
	short int dim_un0;              /* 34 + 2    */
	float pixdim[8];                /* 36 + 32   */
			/*
				pixdim[] specifies the voxel dimensions:
				pixdim[1] - voxel width
				pixdim[2] - voxel height
				pixdim[3] - interslice distance
					..etc
			*/
        float vox_offset;                                   /* 68 + 4    */
	float funused1;                                     /* 72 + 4    */
	float funused2;                                     /* 76 + 4    */
	float funused3;                                     /* 80 + 4    */
        float cal_max;                                      /* 84 + 4    */
        float cal_min;                                      /* 88 + 4    */
	int compressed;                                     /* 92 + 4    */
	int verified;                                       /* 96 + 4    */
	int glmax, glmin;                                   /* 100 + 8   */
	} dime;
	struct data_history
	{
		char descrip[80];
		char aux_file[24];
		char orient;
		char originator[10];
		char generated[10];
		char scannum[10];
		char patient_id[10];
		char exp_date[10];
		char exp_time[10];
		char hist_un0[3];
		int views;
		int vols_added;
		int start_field;
		int field_skip;
		int omax, omin;
		int smax, smin;
	} dh;
};
