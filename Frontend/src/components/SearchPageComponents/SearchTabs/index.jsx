import React from 'react';
import { TabGroup } from 'components/Tab';
import { setSeachCategoryInUrl } from 'shared/utils/common/utils';

const SearchTabs = ({
  value, onChange, searchParams, path,
}) => {
  let categories = ['teams', 'institutions', 'users', 'projects'];
  categories = categories.map((category) => ({
    label: category.toUpperCase(),
    path: `${setSeachCategoryInUrl(path, category)}?${searchParams}`,
  }));
  // Not nice solution but TabGroup works only with integer at the moment.
  const index = categories.findIndex(
    (item) => item.label.toLowerCase() === value,
  );
  return (
    <TabGroup
      value={index}
      setValue={() => onChange()}
      tabsInfo={categories}
    />
  );
};

export default SearchTabs;
