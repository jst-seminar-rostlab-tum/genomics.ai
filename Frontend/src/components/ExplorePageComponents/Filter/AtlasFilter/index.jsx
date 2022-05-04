import React from 'react';
import {
  FormGroup, Stack, Box, Divider,
} from '@mui/material';
import LabeledSelect from 'components/SearchPageComponents/Filter/LabeledSelect';
import FilterSelectTag from '../FilterSelectTag';

const FilterItem = ({ children }) => <Box sx={{ margin: 0.5 }}>{children}</Box>;

// General filter that is needed in all categories
const AtlasFilter = ({ sortBy, onChange, compatibleModels }) => {
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
      <Divider sx={{ margin: 1, width: '100%' }} />
      <FilterItem>
        <FilterSelectTag label="Compatible Models" categories={compatibleModels} onChange={(value) => onChange('compatibleModels', value)} />
      </FilterItem>
    </FormGroup>
  );
};

export default AtlasFilter;
