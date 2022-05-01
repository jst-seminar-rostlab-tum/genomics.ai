import React from 'react';
import { FormGroup } from '@mui/material';
import { Route } from 'react-router-dom';
import { setSeachCategoryInUrl } from 'shared/utils/common/utils';
import { GeneralCard } from 'components/Cards/GeneralCard';
import GeneralFilter from './GeneralFilter';
import TeamsFilter from './TeamsFilter';
// eslint-disable-next-line import/no-extraneous-dependencies
import { Box } from '@mui/system';

const FilterItem = ({ children }) => <Box sx={{ margin: 0.5 }}>{children}</Box>;

// Component storing the necessary filter
function Filter({ path, searchParams, updateQueryParams }) {
  return (
    <GeneralCard>
      <FormGroup>
        <FilterItem>
          <GeneralFilter
            sortBy={searchParams.get('sortBy')}
            onChange={(param, value) => updateQueryParams(param, value)}
          />
        </FilterItem>
        <Route path={setSeachCategoryInUrl(path, 'teams')}>
          <FilterItem>
            <TeamsFilter
              visibility={searchParams.get('visibility')}
              onChange={(param, value) => updateQueryParams(param, value)}
            />
          </FilterItem>
        </Route>
      </FormGroup>
    </GeneralCard>
  );
}

export default Filter;
