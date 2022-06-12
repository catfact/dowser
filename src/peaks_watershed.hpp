#ifndef DOWSER_PEAKS_HPP
#define DOWSER_PEAKS_HPP

#include <array>
#include <limits>
#include <vector>

namespace dowser
{
    namespace peaks
    {

        template <typename T, int N>
        class Watershed
        {
        public:
            class Peak
            {
            public:
                int index;
                int left;
                int right;
                int born;
                int died;
                float persistence;
                T value;

            public:
                Peak(int aIndex) : index(aIndex)
                {
                    born = left = right = index;
                    died = -1;
                }
                float getPersistence(const T *series)
                {
                    if (died >= 0)
                    {
                        return series[born] - series[died];
                    }
                    else
                    {
                        return std::numeric_limits<float>::infinity();
                    }
                }
                void updatePersistence(const T *series)
                {
                    this->persistence = this->getPersistence(series);
                }
                void setValue(T aValue) { value = aValue; }
            };

            static std::vector<Peak> findPeaks(T *series)
            {
                std::vector<int> indices;
                indices.reserve(N);
                for (int i = 0; i < N; ++i)
                {
                    indices.push_back(i);
                }
                std::sort(indices.begin(), indices.end(), [series](int a, int b)
                          { return series[a] > series[b]; });

                std::vector<int> indexToPeak;
                indexToPeak.resize(N);
                std::fill(indexToPeak.begin(), indexToPeak.end(), -1);

                std::vector<Peak> peaks;

                for (int i = 0; i < N; ++i)
                {
                    int idx = indices[i];
                    bool done_l = idx > 0 && indexToPeak[idx - 1] >= 0;
                    bool done_r = idx < (N - 1) && indexToPeak[idx + 1] >= 0;
                    int idx_l = done_l ? indexToPeak[idx - 1] : -1;
                    int idx_r = done_r ? indexToPeak[idx + 1] : -1;

                    // a new peak is born
                    if (!done_l && !done_r)
                    {
                        peaks.push_back(Peak(idx));
                        indexToPeak[idx] = peaks.size() - 1;
                    }

                    // merge to left
                    if (done_l && !done_r)
                    {
                        peaks[idx_l].right += 1;
                        indexToPeak[idx] = idx_l;
                    }

                    // merge to right
                    if (done_r && !done_l)
                    {
                        peaks[idx_r].left -= 1;
                        indexToPeak[idx] = idx_r;
                    }

                    // merge left and right
                    if (done_l && done_r)
                    {
                        if (series[peaks[idx_l].born] > series[peaks[idx_r].born])
                        {
                            peaks[idx_r].died = idx;
                            peaks[idx_l].right = peaks[idx_r].right;
                            indexToPeak[peaks[idx_l].right] = indexToPeak[idx] = idx_l;
                        }
                        else
                        {
                            peaks[idx_l].died = idx;
                            peaks[idx_r].left = peaks[idx_l].left;
                            indexToPeak[peaks[idx_r].left] = indexToPeak[idx] = idx_r;
                        }
                    }
                }
                for (auto &p : peaks)
                {
                    p.updatePersistence(series);
                    p.setValue(series[p.index]);
                }
                std::sort(peaks.begin(), peaks.end(), [](Peak a, Peak b) {
                    return a.persistence > b.persistence;
                });
                return peaks;
            }
        };
    }
}

#endif