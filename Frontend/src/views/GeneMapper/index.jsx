import { Box, Typography, Button } from '@mui/material';
import React, { useCallback, useState } from 'react';
import {
  Redirect, Route, Switch, useRouteMatch,
} from 'react-router-dom';
import GeneMapperState from './GeneMapperState';
import GeneMapperHome from './Home';
import GeneMapperResultView from './Result';
import HeaderView from 'components/general/HeaderView';

function GeneMapper({ sidebarShown }) {
  const paddingL = useCallback(() => (sidebarShown ? '100px' : '350px'), [sidebarShown]);

  const { path, url } = useRouteMatch();

  const headerStyle = {
    color: "#003560",
    fontSize: "1.6rem",
    fontWeight: "bold",
    paddingTop: "20px"
  } 
  const [createOpen, setCreateOpen] = useState(false);

  return (
    <div>
      <HeaderView
        sidebarShown={sidebarShown}
        title="Gene Mapper"
      >
        <Switch>
          <Route exact path={`${path}`}>
            <GeneMapperHome basePath={path} />
          </Route>

          <Route path={`${path}/create`}>
            <GeneMapperState basePath={`${path}/create`} path={`${path}`} />
          </Route>

          <Route path={`${path}/result/:projectId`}>
            <GeneMapperResultView />
          </Route>

          <Route path={`${path}/*`}>
            <Redirect to={`${url}`} />
          </Route>
        </Switch>
      </HeaderView>
    </div>
  );
}

export default GeneMapper;
