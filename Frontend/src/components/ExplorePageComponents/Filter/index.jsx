import { FormGroup, Stack, Box } from '@mui/material';
import React from 'react';
import { Route } from 'react-router-dom';
import { setSeachCategoryInUrl } from 'shared/utils/common/utils';
import { GeneralCard } from 'components/Cards/GeneralCard';
import GeneralFilter from 'components/SearchPageComponents/Filter/GeneralFilter';
import AtlasFilter from './AtlasFilter';

const FilterItem = ({ children }) => <Box sx={{ margin: 0.5 }}>{children}</Box>;

// Component storing the necessary filter
const Filter = ({ path, searchParams, updateQueryParams }) => (
  <GeneralCard>
    <FormGroup>
      <Route path="/explore/atlases">
        <AtlasFilter
          sortBy={searchParams.get('sortBy')}
          onChange={(param, value) => updateQueryParams(param, value)}
        />
      </Route>
      <Route path="/explore/models">
        <FilterItem>
          <GeneralFilter
            sortBy={searchParams.get('sortBy')}
            onChange={(param, value) => updateQueryParams(param, value)}
          />
        </FilterItem>
      </Route>
      <Route path="/explore/datasets">
        <FilterItem>
          <GeneralFilter
            sortBy={searchParams.get('sortBy')}
            onChange={(param, value) => updateQueryParams(param, value)}
          />
        </FilterItem>
      </Route>
    </FormGroup>
  </GeneralCard>
);

export default Filter;
