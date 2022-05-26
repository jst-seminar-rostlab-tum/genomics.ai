import React from 'react';
import { FormGroup, Box } from '@mui/material';
import { Route } from 'react-router-dom';
import { setSeachCategoryInUrl } from 'shared/utils/common/utils';
import { GeneralCard } from 'components/Cards/GeneralCard';
import GeneralFilter from './GeneralFilter';
import TeamsFilter from './TeamsFilter';
import AtlasFilter from 'components/ExplorePageComponents/Filter/AtlasFilter';

const FilterItem = ({ children }) => <Box sx={{ margin: 0.5 }}>{children}</Box>;

// Component storing the necessary filter
function Filter({ path, searchParams, updateQueryParams }) {
  return (
    <GeneralCard>
      <FormGroup>
        <Route path={[setSeachCategoryInUrl(path, 'teams'),
          setSeachCategoryInUrl(path, 'institutions'),
          setSeachCategoryInUrl(path, 'models')]}
        >
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
        </Route>
        <Route path={setSeachCategoryInUrl(path, 'atlases')}>
          <AtlasFilter
            sortBy={searchParams.get('sortBy')}
            onChange={(param, value) => updateQueryParams(param, value)}
          />
        </Route>
        <Route path={setSeachCategoryInUrl(path, 'users')}>
          <FilterItem>
            <GeneralFilter
              sortBy={searchParams.get('sortBy')}
              onChange={(param, value) => updateQueryParams(param, value)}
              sortItems={[{ label: 'Name', value: 'name' }]}
              defaultValue="name"
            />
          </FilterItem>
        </Route>
      </FormGroup>
    </GeneralCard>
  );
}

export default Filter;
