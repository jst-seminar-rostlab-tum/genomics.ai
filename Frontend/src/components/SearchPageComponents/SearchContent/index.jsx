import React from 'react';

import {
  useRouteMatch, Switch, Route,
} from 'react-router-dom';
import SearchResultList from '../SearchResultList';
import TeamCard from '../SearchResultList/TeamCard';
import InstitutionCard from '../SearchResultList/InstitutionCard';
import UserCard from '../SearchResultList/UserCard';
import ResultStatus from '../ResultStatus';
import { setSeachCategoryInUrl } from 'shared/utils/common/utils';
import AtlasesGrid from 'components/Grids/AtlasesGrid';
import ModelsGrid from 'components/Grids/ModelsGrid';
import ExploreRoutes from 'components/ExplorePageComponents/ExploreRoutes';

// wrapper component to display the searched items
function SearchContent({
  searchResult, searchCategory, searchedKeyword, user, fetchSearchHandler,
}) {
  const { path } = useRouteMatch();

  const renderSearchResultsList = (listItemWrapper) => (
    <SearchResultList
      listItemWrapper={listItemWrapper}
      searchResult={searchResult}
      user={user}
      fetchSearchHandler={fetchSearchHandler}
    />
  );

  const atlases = <AtlasesGrid atlases={searchResult} searchedKeyword={searchedKeyword} path="/sequencer/search" isSearchPage/>;
  const models = <ModelsGrid models={searchResult} searchedKeyword={searchedKeyword} path="/sequencer/search" isSearchPage/>;
  return (
    <>
      <ResultStatus
        count={searchResult.length}
        searchedEntity={searchCategory}
        searchedKeyword={searchedKeyword}
      />
      <Switch>
        <Route path={setSeachCategoryInUrl(path, 'teams')}>
          {renderSearchResultsList(TeamCard)}
        </Route>
        <Route path={setSeachCategoryInUrl(path, 'institutions')}>
          {renderSearchResultsList(InstitutionCard)}
        </Route>
        <Route path={setSeachCategoryInUrl(path, 'users')}>
          {renderSearchResultsList(UserCard)}
        </Route>
        <Route>
          <ExploreRoutes atlases={atlases} models={models} path="/sequencer/search" />
        </Route>
      </Switch>
    </>
  );
}

export default SearchContent;
