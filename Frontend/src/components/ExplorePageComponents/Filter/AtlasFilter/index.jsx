import React from 'react';
import {
  FormGroup, Box,
} from '@mui/material';
import LabeledSelect from 'components/SearchPageComponents/Filter/LabeledSelect';

const FilterItem = ({ children }) => <Box sx={{ margin: 0.5 }}>{children}</Box>;

// General filter that is needed in all categories
const AtlasFilter = ({ sortBy, onChange }) => {
  const sortItems = [
    { label: 'Name', value: 'name' },
    { label: 'Last updated', value: 'lastUpdated' },
    { label: 'Number of Cells', value: 'numberOfCells' },
  ];

  return (
    <FormGroup sx={{ alignItems: 'start' }}>
      <FilterItem>
        <LabeledSelect
          label="Sort by"
          value={sortBy}
          defaultValue="name"
          onChange={(event) => onChange('sortBy', event.target.value)}
          items={sortItems}
        />
      </FilterItem>
    </FormGroup>
  );
};

export default AtlasFilter;
