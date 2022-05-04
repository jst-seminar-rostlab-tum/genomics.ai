import { FormGroup, Stack, Box } from '@mui/material';
import React from 'react';
import { Route } from 'react-router-dom';
import { setSeachCategoryInUrl } from 'shared/utils/common/utils';
import { GeneralCard } from 'components/Cards/GeneralCard';
import GeneralFilter from 'components/SearchPageComponents/Filter/GeneralFilter';

const FilterItem = ({ children }) => <Box sx={{ margin: 0.5 }}>{children}</Box>;

// Component storing the necessary filter
const Filter = ({ searchParams, updateQueryParams }) => (
  <GeneralCard>
    <FormGroup>
      <FilterItem>
        <GeneralFilter
          sortBy={searchParams.get('sortBy')}
          onChange={(param, value) => updateQueryParams(param, value)}
        />
      </FilterItem>
    </FormGroup>
  </GeneralCard>
);

export default Filter;
