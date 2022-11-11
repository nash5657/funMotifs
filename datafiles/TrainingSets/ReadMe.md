This directory must have a specified structure due to the function get_trainings_data_dirs
and get_param_weight in WeightFeatures.py \
This simplifies the parsing of the data for the user.\
The structure is as follows:\

>TrainingSets
>>tile_prom_region
> >>*cell_name1*
> >>>input_data
> >>>>infile.txt
> >>>
> >>>output_data
> >>
> >prom_unactive
> >>*cell_name2*
> >>>input_data
> >>>>*cell_name2*.bed
> >>>
> >>>output_data
> >>
> >other_active_region
> >>*cell_name3*
> >>>input_data
> >>>>infile_coordinates.bed\
> >>>>infile_expression.txt
> >>>
> >>>output_data

Data in other format will be ignored by the code.