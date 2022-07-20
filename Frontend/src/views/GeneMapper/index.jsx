import { Container, Box, Typography } from '@mui/material';
import React, { useCallback, useState } from 'react';
import {
  Redirect, Route, Switch, useRouteMatch,
} from 'react-router-dom';
import GeneMapperState from './GeneMapperState';
import NavBar from 'components/NavBar';
import GeneMapperHome from './Home';
import GeneMapperResultView from './Result';
import HeaderView from 'components/general/HeaderView';
import HelpIcon from '@mui/icons-material/Help';
import { Modal } from 'components/Modal';
import { colors } from 'shared/theme/colors';

/*
params: sidebarShown, loggedIn
Logged-in version of the GeneMapper. loggedIn is set to true by default
*/
function GeneMapper({ sidebarShown, loggedIn = true }) {
  const { path, url } = useRouteMatch();

  const headerStyle = {
    color: '#003560',
    fontSize: '1.6rem',
    fontWeight: 'bold',
    paddingTop: '20px',
  };
  const [showScarchesInfo, setShowScarchesInfo] = useState(false);

  return (
    <div>
      <HeaderView
        sidebarShown={sidebarShown}
        title="ScArches Gene Mapper"
        rightOfTitle={(
          <Box sx={{ display: 'flex', alignItems: 'center' }}>
            &nbsp;
            <HelpIcon
              onClick={() => setShowScarchesInfo(true)}
              color="primary"
              sx={{
                cursor: 'pointer',
              }}
            />
          </Box>
        )}
        loggedIn={loggedIn}
      >
        <Switch>
          <Route exact path={`${path}`}>
            {/* should mostly be fine, todo: make sure the paths work properly */}
            <GeneMapperHome basePath={path} loggedIn={loggedIn} />
          </Route>

          <Route path={`${path}/create`}>
            {/* should mostly be fine, todo: make sure the paths work properly */}
            <GeneMapperState path={`${path}`} loggedIn={loggedIn} />
          </Route>

          <Route path={`${path}/result/:projectId`}>
            <GeneMapperResultView loggedIn={loggedIn} />
          </Route>

          <Route path={`${path}/*`}>
            <Redirect to={`${url}`} />
          </Route>
        </Switch>
      </HeaderView>
      <Modal
        isOpen={showScarchesInfo}
        setOpen={setShowScarchesInfo}
        children={(
          <Container>
            <Box sx={{
              display: 'flex',
              flexDirection: 'column',
              alignItems: 'flex-start',
              justifyContent: 'space-between',
            }}
            >
              <Typography sx={{ fontSize: '36px', fontWeigth: 700 }}>What is ScArches?</Typography>
              <Typography sx={{ fontSize: '16px', fontWeight: 300, paddingBottom: '16px' }}>
                <p>
                  ScArches is a novel deep learning model that enables mapping query to reference datasets. The model allows the user to construct single or multi-modal (CITE-seq) references as well as classifying unlabelled query cells.
                  On ArchMap, currently only mapping of query to reference datasets is offered.

                </p>
                <p>
                  For further information, please check out the ScArches paper:
                  <a href="https://www.nature.com/articles/s41587-021-01001-7" target="_blank" rel="noreferrer">https://www.nature.com/articles/s41587-021-01001-7</a>
                </p>
              </Typography>
            </Box>
          </Container>
        )}
      />
    </div>
  );
}

export default GeneMapper;
