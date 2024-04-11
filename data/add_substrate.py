import pandas as pd

if __name__ == "__main__":
    d1 = pd.read_csv('growth_rate_ds.csv', index_col=0)
    d2 = pd.read_csv('extracted_growth_rates.csv', index_col=0)

    d1 = d1[['image_id', 'medium']].drop_duplicates(ignore_index=True)
    print(d1)

    # print(d2)

    image_ids = d1['image_id']
    media = d1['medium']

    _map = {}
    for i,image_id in enumerate(image_ids):
        _map[image_id] = media[i]

    _map[26185] = 'PCA-Gluc'
    _map[26186] = 'PCA-Gluc'
    _map[26187] = 'PCA-Gluc'
    _map[26188] = 'PCA-Gluc'
    _map[26189] = 'PCA-Gluc'
    _map[26190] = 'PCA-Gluc'
    _map[26191] = 'PCA-Gluc'
    _map[26192] = 'PCA-Gluc'
    _map[26193] = 'PCA-Gluc'
    _map[26194] = 'PCA-Gluc'
    _map[26195] = 'PCA-Gluc'

    print(_map)

    d2['medium'] = d2['image_id']
    print(d2['image_id'])
    d2['medium'] = d2['medium'].apply(lambda x: _map[x] if x in _map else 0)
    print(d2)
    d2.to_csv('output.csv')
    # res = d2.join(d1, how="outer", on="image_id")
    #
    # print(res)
