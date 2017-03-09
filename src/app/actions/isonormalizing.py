from statistics import median, StatisticsError


def get_medians(channels, ratios):
    ch_medians = {}
    for ix, channel in enumerate(channels):
        try:
            ch_medians[channel] = median([x[ix] for x in ratios
                                          if x[ix] != 'NA'])
        except StatisticsError:
            # channel is empty, common in protein quant but not in normalizing
            ch_medians[channel] = 'NA'
    return ch_medians
