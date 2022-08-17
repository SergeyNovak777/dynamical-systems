# dynamical-systems

## Containt repository:
 * Hindmarrh Rose (папка для системы Hindmarrh Rose)
   * chaos_EV_map
     * EVA_spikes.npy (массив ЭС. std, mean принимают спайки x_sum)
     * EVA_without_spikes.npy (массив ЭС. std, mean принимают x_sum)
     * LE.npy (Массив старшего ляпуновского показателя)
     * Hindmarsh–Rose Map Lyapunov.ipynb (для карт ЭС и старшего ляпуновского показателя)
     
   * attractors 
     * Hindmarsh Rose poincaremap.ipynb (для отображения Пуанкаре)
     * FixedPoints solo element.ipynb (Неподвижные точки для одиночного элемента)
     * FixedPoints.ipynb (Неподвижные точки для двух связанных HR)
     * Modes.ipynb (Режимы)
     * basins code from doc ( бассейны аттракторов, примеры из документации)
     * basins of attractors chaos.ipynb ( бассейны аттракторов в хаосе)
     * basins of attractors quasi.ipynb ( бассейны аттракторов при k = 0.0)
     * mode comparison.ipynb ( сравнение квазипериодического и хаотического режима)
     
   * probablity density function
     * result ( содержит папки в которых сохранены изображения для PDF, timeseries, lyapunov spectrum)
       * small_step_0.001
         * dots
         * full
         * pdf
         * timeseries
         * condition_HR_array.jld (массив из начальных условий. Julia)
         * condition_HR_array.npy (массив из начальных условий. Python)
         * spectrum_HR_array.jld (массив из значений спектра Ляпунова. Julia)
         * spectrum_HR_array.npy (массив из значений спектра Ляпунова. Python)
         
     * Hindmarsh Rose PDF.ipynb ( вычисление PDF, lyapunov specrtrum с их сохранением)
     * Benchmark.ipynb ( Оптимизация системы, замеры времени)
     * PDF three value k.ipynb ( детальное рассмотрение 3 значений k из папки small_step_0.001)
  * Rulkov (Папка для системы Рулькова)
    
#### Система в большинстве блокнотов не оптимизирована. Версии оптимизированной системы находятся в "dynamical-systems/Hindmarrh Rose/probablity density function/Benchmark.ipynb".
#### Рекомендуется использовать вторую версию оптимизации вместе с SA. Для этого необходимо установить пакет StaticArrays.
 
