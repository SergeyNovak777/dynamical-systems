# dynamical-systems

## Containt repository:
 * Hindmarrh Rose (папка для системы Hindmarrh Rose)
   * chaos_EV_map 
     * Hindmarsh–Rose Map Lyapunov.ipynb (для карт ЭС и старшего ляпуновского показателя)
   * pouncare map 
     * Hindmarsh Rose poincaremap.ipynb (для отображения Пуанкаре)
   * probablity density function
     * result ( содержит папки в которых сохранены изображения для PDF, timeseries, lyapunov spectrum)
       * small_step_0.001
         * dots
         * full
         * pdf
         * timeseries
         * condition_HR_array.jld (массив из начальных условий. Сохранен в формате для Julia)
         * condition_HR_array.npy (массив из начальных условий. Сохранен в формате для Python)
         * spectrum_HR_array.jld (массив из значений спектра Ляпунова. Сохранен в формате для Julia)
         * spectrum_HR_array.npy (массив из значений спектра Ляпунова. Сохранен в формате для Julia)
         
     * FixedPoints.ipynb (Поиск неподвижных точек для исходной системы)
     * Hindmarsh Rose Benchmark.ipynb ( Оптимизация системы, замеры времени)
     * Hindmarsh Rose PDF three value k.ipynb ( детальное рассмотрение 3 значений k из папки HR_small_saved)
     * Hindmarsh Rose PDF.ipynb ( вычисление PDF, lyapunov specrtrum с их сохранением)
    
