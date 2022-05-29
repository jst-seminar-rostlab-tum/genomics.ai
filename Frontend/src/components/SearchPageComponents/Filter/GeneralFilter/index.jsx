import React from 'react';
import LabeledSelect from '../LabeledSelect';

// General filter that is needed in all categories
const GeneralFilter = ({ sortBy, onChange }) => {
  const defaultValue = 'name';

  const sortItems = [
    { label: 'Name', value: 'name' },
    { label: 'Last updated', value: 'updatedAt' },
  ];

  const sortItem = sortItems.find((item) => item.value === sortBy);

  return (
    <LabeledSelect
      label="Sort by"
      value={sortItem && sortItem.value}
      defaultValue={defaultValue}
      onChange={(event) => onChange('sortBy', event.target.value)}
      items={sortItems}
    />
  );
};

export default GeneralFilter;
