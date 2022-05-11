import React from 'react';

import { useRouteMatch, Switch, Route } from 'react-router-dom';
import SearchResultList from '../SearchResultList';
import TeamCard from '../SearchResultList/TeamCard';
import InstitutionCard from '../SearchResultList/InstitutionCard';
import UserCard from '../SearchResultList/UserCard';
import ResultStatus from '../ResultStatus';
import { setSeachCategoryInUrl } from 'shared/utils/common/utils';

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
      </Switch>
    </>
  );
}

export default SearchContent;
