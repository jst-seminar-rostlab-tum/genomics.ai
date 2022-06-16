import React from 'react';
import { TabGroup } from 'components/Tab';
import { setSeachCategoryInUrl } from 'shared/utils/common/utils';

const SearchTabs = ({
  value, searchParams, path, onChange = () => {},
}) => {
  let categories = ['atlases', 'models', 'teams', 'institutions', 'users'];
  const keyword = searchParams.get('keyword');
  const newParams = keyword ? `?keyword=${keyword}` : '';
  categories = categories.map((category) => ({
    label: category.toUpperCase(),
    path: `${setSeachCategoryInUrl(path, category)}${newParams}`,
  }));
  // Not nice solution but TabGroup works only with integer at the moment.
  const index = categories.findIndex(
    (item) => item.label.toLowerCase() === value,
  );
  return (
    <TabGroup
      value={index}
      onValueChange={() => onChange()}
      setValue={() => onChange()}
      tabsInfo={categories}
    />
  );
};

export default SearchTabs;
